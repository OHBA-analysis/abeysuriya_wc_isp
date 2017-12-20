function [sol,sys] = integrate_wc( unit, coupl, delay, input, hist, tspan, dt )
    % Integrate network of Wilson-Cowan neural masses
    %
    % INPUTS:
    %
    % - unit : Node parameters. For a network with N regions, there are 2N
    %   nodes. This is a struct where each field is a vector of length 2N. The
    %   fields required are:
    %       - mu : (sigmoid mean)
    %       - sigma : (sigmoid width, >0)
    %       - qmax : (maximum firing rate, >0)
    %       - tau : (time constant, s^-1, >0)
    %       - ce,ci : local excitatory (>=0) / inhibitory (<=0) coupling Note
    %         that cee <= unit.ce(2p), cei <= unit.ce(2p+1), cie <=
    %         unit.ci(2p), cii <= unit.ci(2p+1)
    %       - r : refractory "period" (>=0)
    %       - isp_rate : On inhibitory nodes (i.e. odd numbered entries), this
    %         is 1/tau_isp
    %       - isp_target : On inhibitory notes, this is rho
    % - coupl :  NxN matrix of couplings (>=0). This matrix should have zeros
    %   on the diagonal (i.e. delayed self-connections not supported) and all
    %   long-range couplings (cf matrix coupl) should be non-negative.
    % - delay : NxN matrix of delays (>=0)
    % - input :  Node inputs. Struct array with 2N entries, with fields
    %   {time,vals}, where:
    %       - time : 1xNt vector with time values
    %       - vals : 1xNt vector with input values The inputs are interpolated
    %         onto the integration grid (half-steps for RK4). Nodes thus do
    %         not all require the same input times (e.g. constant inputs can
    %         be represented with 3-element vectors)
    % -  hist : Node history. Similarly to inputs, structure with fields
    %    {time,vals}, where vals is a 2NxNh matrix. The history span should be
    %    long enough to accomodate the largest delay.
    % - tspan : Time-span of simulation (starts at hist.time(end))
    % - dt : Integration time-step. The smallest positive delay should be larger than the time-step.
    %
    % OUTPUTS:
    %
    % - sol : Struct with fields {time,vals}, where vals is a 2NxNt matrix.
    % - sys : Struct with fields {unit,edge,input} = system used for integration.
    %
    % Romesh Abeysuriya Nov17
    % Jonathan Hadida Nov17 @fmrib @ohba

    [Nu,Nn] = check_inputs( unit, coupl, delay, input, hist, tspan, dt );
    
    % time-frame
    Tini = hist.time(end);
    Tend = Tini + tspan;
    Tpts = Tini : dt : Tend;
    Nt = numel(Tpts);
    
    % integrate history into solution
    hist = resample_history( hist, dt );
    check_hist_delay( hist, delay, dt );
    
    % allocate solution (including history)
    sol.time = [ hist.time, Tpts(2:end) ];
    sol.vals = [ hist.vals, nan(Nn,Nt-1) ];
    sol.cie = nan(Nn/2,size(sol.vals,2));

    Nh = numel(hist.time); % start index
    
    % vectorise unit fields
    unit = structfun( @(v) v(:), unit, 'UniformOutput', false );
    
    % prepare inputs for integration
    input_interp.time = Tini:dt/2:Tend;
    vals = zeros(Nn,length(input_interp.time));
    for j = 1:length(input)
        vals(j,:) = pchip( input(j).time, input(j).vals , input_interp.time);
    end
    input_interp.vals = vals;
    
    % process couplings and delays
    edge = process_network( coupl, delay );
    [edge,cie_index] = add_local( edge, unit );
    edge = precompute_edges(edge,Nn,dt);

    sol.cie(:,Nh) = edge.instant.c(cie_index);

    % wrap-up system so it can be output with solution
    sys.unit = unit;
    sys.edge = edge;
    sys.input = input;
    
    % iterate on time-points
    for i = 1:Nt-1
        
        % index of current and next timepoint in sol
        nxt = Nh + i; 
        cur = nxt - 1;
        
        % current data
        ti = sol.time(cur);
        xi = sol.vals(:,cur);

        % update cache
        cache = make_cache( edge, sol, cur, dt );
        
        % compute intermediate RK4 points
        k1 = wc_eval( unit, input_interp.vals(:,1+2*(i-1)+0), cache, edge, ti, xi, dt );
        k2 = wc_eval( unit, input_interp.vals(:,1+2*(i-1)+1), cache, edge, ti + dt/2, xi + k1*dt/2, dt );
        k3 = wc_eval( unit, input_interp.vals(:,1+2*(i-1)+1), cache, edge, ti + dt/2, xi + k2*dt/2, dt );
        k4 = wc_eval( unit, input_interp.vals(:,1+2*(i-1)+2), cache, edge, ti + dt, xi + k3*dt, dt );
        
        % final estimate
        yi = xi + (dt/6) * ( k1 + 2*k2 + 2*k3 + k4 );
        
        % save next point
        sol.vals(:,nxt) = yi;
        edge.instant.c(cie_index) = edge.instant.c(cie_index) - dt.*unit.isp_rate(2:2:end).*yi(2:2:end).*(yi(1:2:end)-unit.isp_target(2:2:end));
        sol.cie(:,nxt) = edge.instant.c(cie_index);

    end
    
    % extract only those timepoints after history
    sol.time = sol.time(Nh:end);
    sol.vals = sol.vals(:,Nh:end);
    sol.cie = sol.cie(:,Nh:end);

end

function y = wc_eval( unit, input, cache, edge, t, x, dt )
%
% Evaluate the Wilson-Cowan system
%

    n = numel(x);
    
    % interpolate delayed values and inputs
    w = (cache.tr - t) / dt; %assert(w >= 0);
    y = w * cache.yl + (1-w) * cache.yr + input;
    
    % add local contributions
    y = y + edge.instant.accum_matrix * (edge.instant.c.*x(edge.instant.src));

    % apply sigmoid
    y = (y - unit.mu) ./ unit.sigma;
    y = -x + (1 - unit.r .* x) .* unit.qmax ./ ( 1 + exp(-y) );
    y = y ./ unit.tau;
    
end

function cache = make_cache( edge, sol, k, dt )
%
% Cache delayed terms for faster integration.
%

    n = size(sol.vals,1);
    
    cache.tl = sol.time(k);
    cache.tr = cache.tl+dt;
    
    if isempty(edge.delayed.d)
        cache.yl = zeros(n,1);
        cache.yr = zeros(n,1);
    else
        r = edge.delayed.src + (k-edge.delayed.steps)*n;
        
        s2 = sol.vals(r-2*n);
        s1 = sol.vals(r-n);
        s0 = sol.vals(r);
        
        yl = edge.delayed.step_remainder .* s2  +  (1-edge.delayed.step_remainder) .* s1; % linear interpolation
        yr = edge.delayed.step_remainder .* s1  +  (1-edge.delayed.step_remainder) .* s0; 
            
        cache.yl = edge.delayed.accum_matrix * (yl);
        cache.yr = edge.delayed.accum_matrix * (yr);

    end

end

function [Nu,Nn] = check_inputs( unit, coupl, delay, input, hist, tspan, dt )
%
% Basic checks to make sure everything is in order.
%

    ufields = {'mu','sigma','qmax','tau','ce','ci','r','isp_rate','isp_target'};
    assert( check_struct( unit, ufields ), ...
        'Unit should be a structure with fields: %s', strjoin(ufields,', ') );

    Nn = numel(unit.mu);
    Nu = Nn / 2;
    assert( 2*Nu == Nn, 'Number of nodes should be even.' );
    
    assert( all(structfun( @numel, unit ) == Nn), 'All unit fields should be the same length.' );
    assert( all(unit.ci <= 0), 'Local inhibitory couplings should be non-positive.' );
    assert( all(unit.ce >= 0), 'Local excitatory couplings should be non-negative.' );
    
    assert( check_square_noneg(coupl,Nu), 'Coupling should be a %dx%d matrix with non-negative entries.', Nu, Nu );
    assert( check_square_noneg(delay,Nu), 'Delay should be a %dx%d matrix with non-negative entries.', Nu, Nu );
    
    assert( isstruct(input) && all(isfield(input,{'time','vals'})), ...
        'Input should be a struct with fields {time,vals}' );
    assert( check_tcs(hist, Nn), ...
        'Hist should be a struct with fields {time,vals}, and vals should have %d rows.', Nn );
    
    assert( check_posnum(tspan), 'tspan should be a positive number.' );
    assert( check_posnum(dt), 'dt should be a positive number.' );
    
end

function check_hist_delay( hist, delay, dt )
%
% Make sure that history time-courses are long-enough to accomodate the largest delay,
% and that the smallest positive delay is larger than the time-step.
%
% ISSUE:
%   Actually, the largest delay may be ignored if the associated coupling is null.
%   This is a minor issue which may cause an error even though the integration would be possible.
%

    dmax = max(delay(:));
    dmin = min(delay(delay > 0));
    hspan = hist.time(end) - hist.time(1);
    
    assert( all(dmin >= dt), 'Minimum positive delay (%.2f ms) should be larger than the timestep (%.2f ms).', ...
        1000*dmin, 1000*dt );
    assert( all(dmax <= hspan), 'History is too short (%.2f ms) to accomodate the largest delay (%.2f ms).', ...
        1000*hspan, 1000*dmax );

end

function hist = resample_history( hist, dt )
%
% Make sure the history time-courses are arithmetically sampled with the same time-step as the output.
%
% NOTE:
%   There is NO EXTRAPOLATION that takes place here, so make sure the history is long enough to accomodate
%   the largest delay. If (hist.time(end) - hist.time(1))/dt is not an integer, then the resampled history
%   will be shorter than whatever was input manually.
%
    t = fliplr( hist.time(end) : (-dt) : hist.time(1) );
    hist.time = t;
    hist.vals = pchip(hist.time,hist.vals,t);
    
end

function edge = process_network( coupl, delay )
%
% edge is a struct with fields {delayed,instant}, each with fields:
%   src  source node
%   dst  destination node
%     c  non-zero coupling
%     d  corresponding delay (not for instant)
%
% NOTE:
%   Only positive couplings are retained after this processing.
%   You can change the value of THRESH below to adjust how "positive" the minimum coupling should be.
%   

    THRESH = 1e-12; % couplings smaller than this will be ignored
    
    % force diagonals to 0
    Nu = size(coupl,1);
    dg = 1:(Nu+1):(Nu*Nu); % diagonal indices
    coupl(dg) = 0;
    delay(dg) = 0;
    
    % extract edges with non-zero coupling
    [src,dst] = find(coupl > THRESH);
    k = sub2ind( size(coupl), src, dst );
    c = coupl(k);
    d = delay(k);
    z = d > eps;
    
    % long-range couplings are E-E only
    src = 2*src-1;
    dst = 2*dst-1;
    
    % instantaneous edges
    edge.instant.c   = c(~z);
    edge.instant.src = src(~z);
    edge.instant.dst = dst(~z);
    
    % delayed edges
    edge.delayed.c = c(z);
    edge.delayed.d = d(z);
    edge.delayed.src = src(z);
    edge.delayed.dst = dst(z);
    
end

function [edge,cie_index] = add_local( edge, unit )
%
% Add edges corresponding to local couplings.
% These edges have 0 delay, and target either themselves, or the other local population.
%
%    unit.ce        unit.ci
%   1->Cee->1      2->Cie->1
%   1->Cei->2      2->Cii->2
%   3->Cee->3      4->Cie->3
%   3->Cei->4      4->Cii->4
%   5->Cee->5      6->Cie->5
%  src ... dst
%
% cie_index returns the indices of cie in the instant edge list

    n = numel(unit.ce);
    e = [1;1] * (1:2:n); % double each index
    e = e(:);
    i = e+1;
    k = 1:n;
    k = k(:);

    edge.instant.c   = [ edge.instant.c; unit.ce; unit.ci ];
    edge.instant.src = [ edge.instant.src; e; i ];
    edge.instant.dst = [ edge.instant.dst; k; k ];

    cie_index = size(edge.instant.c,1)-size(unit.ci,1)+1:2:size(edge.instant.c,1);

end

function edge = precompute_edges(edge,Nn,dt)
    
    if ~isempty(edge.delayed.d)
        edge.delayed.steps = floor(edge.delayed.d / dt); % Number of whole steps back
        edge.delayed.step_remainder = mod(edge.delayed.d / dt,1); % Amount by which each delayed edge needs to be interpolated relative to a time step
        
        x = zeros(Nn,length(edge.delayed.dst));
        for j = 1:Nn
            x(j,find(edge.delayed.dst==j))=1;
        end
        x = bsxfun(@times,x,edge.delayed.c.');
        edge.delayed.accum_matrix = sparse(x);
    end
    
    x = zeros(Nn,length(edge.instant.dst));
    for j = 1:Nn
        x(j,find(edge.instant.dst==j))=1;
    end
    % Note - cannot premultiply local weights here because they might change due to plasticity
    edge.instant.accum_matrix = sparse(x);
end

% Simplified checks for basic types
function ok = check_posnum(x)
    ok = isscalar(x) && isreal(x) && (x > eps);
end
function ok = check_square_noneg(x,n)
    ok = ismatrix(x) && isreal(x) && all(size(x) == n) && all(x(:) >= 0);
end
function ok = check_struct(x,f)
    ok = isscalar(x) && isstruct(x) && all(isfield(x,f));
end
function ok = check_tcs(x,nrows)
    ok = check_struct(x,{'time','vals'}) && (size(x.vals,1) == nrows);
end