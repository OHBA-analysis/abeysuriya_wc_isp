classdef WilsonCowan < ra.model.base
    
    properties 
        nodes_per_unit = 2;
        states_per_node = 1;
        mx_helper = 'wc_network_rk2'
    end

    methods
        
        function self = WilsonCowan(unit_params,net)
            if nargin == 0 % Default constructor used for preallocation
                return
            end

            self.param = struct();
            self.options = struct();
            self.version = 1;
            self.network = net; 

            if length(unit_params) == 1
                self.units = repmat(unit_params,net.n_nodes,1);
            else
                assert(length(unit_params)==net.n_nodes)
                self.units = unit_params;
            end
            
            % Default params
            self.param.velocity              = 1;
            self.param.coupling              = [1,0]; % Long range EE and EI coupling scaling
            self.param.noise_seed            = 1; % Seed to be used when generating noise timecourses
            self.param.init                  = []; % Initial conditions, must be vector of length self.n_states if specified
            self.param.init_seed             = 2; % Seed to be used when generating random initial conditions
            self.param.round_delays          = false; % Should delays be rounded? This probably matters when noise is present
            self.param.enable_isp            = 0; % Ignored by wc_plastic, plasticity calculation is always carried out (but negligible cost if rate is zero)

            % Default options
            self.options.tstep               = 5e-5;
            self.options.tspan               = 30;
            self.options.cache_delays        = false; % approximate delays at each step
            self.options.induction           = false; 
            self.options.use_pchip           = false;  
            self.options.fixed_step_size     = true; % used to optimise computation of delayed terms
            self.options.compression_epsilon = 1e-6; % coupling values with lower magnitude will be ignored
            self.options.num_threads         = 1;
            self.options.noise_outside       = false;
            self.options.record_feedback     = false;
            self.options.verbose             = false;

        end
        
        % Return common result subset variables
        function ts = excitatory(self), ts = self.result.select_signals(1:self.nodes_per_unit*self.states_per_node:self.n_states); end % First output variable for each unit
        function ts = inhibitory(self), ts = self.result.select_signals(2:self.nodes_per_unit*self.states_per_node:self.n_states); end % Second output variable for each unit
        function ts = hip(self), ts = self.nsys_output.a_hip; end
        function ts = isp(self), ts = self.nsys_output.a_isp; end

        % Return common analysis measures
        function v = cie(self); v = arrayfun(@(x) x.cie,self.units); end % Array of cie values from the units
        function v = iwe(self); v = sum(self.excitatory.vals.*self.inhibitory.vals)./sum(self.inhibitory.vals); end % I-weighted mean excitatory activity
        function d = mean_delay(self); [~,d] = self.network.netmats; d = mean(d(logical(triu(ones(size(d)),1))))./self.param.velocity; end

        function plot(self)
            plot(self.excitatory);
            ylabel('Excitatory activity')
        end

        function h = analyse(self)
            % Put a simple function here to perform an alternate analysis to 'plot'
            % This function is expected to change over time to whatever is convenient
            %ra.analysis.main([3 0.25 0],75).analyse(self);
            % Content below maps to self.plot() for backwards compatibility
            h = figure;
            userdata.tg = uitabgroup(h, 'Position', [0 0 1 1]);
            set(h,'UserData',userdata,'Position',[440   291   954   507]);
            a = uitab(userdata.tg,'Title','Raw TS');
            b = uicontainer(a);
            ax = axes(b);
            plot(self);
        end

        function run(self,varargin)
            % WilsonCowan model also handles the extra NSYS outputs for HIP and ISP
            arg = inputParser;
            arg.addParameter('solver','matlab');
            arg.KeepUnmatched = true;
            arg.parse(varargin{:});

            opt        = ra.utils.kwArgs(varargin{:});
            burn       = opt.get('burn',[0 0]); % Two elements corresponding to seconds to burn from the start and end
            downsample = opt.get('downsample',false); % downsample to a given sampling rate 

            switch arg.Results.solver
                case 'nsys' 
                    fprintf('Using NSYS solver\n')
                    run@ra.model.base(self,varargin); % Perform a standard simulation run
                    self.nsys_output.a_hip = ra.TimeSeries(self.nsys_output.a_hip.time,self.nsys_output.a_hip.vals);
                    self.nsys_output.a_isp = ra.TimeSeries(self.nsys_output.a_isp.time,self.nsys_output.a_isp.vals);
                case 'matlab'
                    fprintf('Using Matlab solver\n')
                    if any(arrayfun(@(x) x.E.hip_rate,self.units))
                        error('HIP not supported by Matlab solver');
                    end
                    
                    unit = struct;
                    m = self.get_mex_inputs(self.options.tspan + sum(burn));
                    unit.mu = arrayfun(@(x) x.S.mu,m.system.nodes);
                    unit.sigma = arrayfun(@(x) x.S.sigma,m.system.nodes);
                    if isfield(m.system.nodes(1).S,'qmax')
                        unit.qmax = arrayfun(@(x) x.S.qmax,m.system.nodes);
                    else
                        unit.qmax = ones(size(unit.mu));
                    end
                    unit.tau = arrayfun(@(x) x.tau,m.system.nodes);
                    unit.r = arrayfun(@(x) x.rp,m.system.nodes);
                    unit.isp_rate = arrayfun(@(x) x.isp_rate,m.system.nodes);
                    unit.isp_target = arrayfun(@(x) x.isp_target,m.system.nodes);

                    unit.ce(1:2:self.n_nodes) = arrayfun(@(x) x.cee,self.units);
                    unit.ce(2:2:self.n_nodes) = arrayfun(@(x) x.cei,self.units);
                    unit.ci(1:2:self.n_nodes) = arrayfun(@(x) x.cie,self.units);
                    unit.ci(2:2:self.n_nodes) = arrayfun(@(x) x.cii,self.units);

                    [cv,dv] = self.network.netmats;
                    coupl = cv * self.param.coupling(1);
                    if self.param.coupling(2) ~= 0
                        error('Matlab integrator does not support long range E->I connections');
                    end
                    delay = dv./self.param.velocity;

                    input = struct('time',arrayfun(@(x) x.P.time,m.system.nodes,'UniformOutput',false),'vals',arrayfun(@(x) x.P.vals,m.system.nodes,'UniformOutput',false));

                    hist.time = m.problem.history.time.';
                    hist.vals = m.problem.history.vals.';
                    tspan = self.options.tspan+sum(burn);
                    dt = self.options.tstep;
                    
                    sol = integrate_wc( unit, coupl, delay, input, hist, tspan, dt );
                    
                    self.result = ra.TimeSeries(sol.time,sol.vals);

                    self.nsys_output.a_hip = ra.TimeSeries(sol.time,bsxfun(@times,arrayfun(@(x) x.E.S.mu,self.units),ones(1,length(sol.time))));
                    self.nsys_output.a_isp = ra.TimeSeries(sol.time,sol.cie);
                    
                    % Replicate downsampling and burn from model.base()
                    if downsample
                        self.result.resample(downsample);
                    end

                    if any(burn)
                        self.result.burn(burn(1));
                        self.result.rem_after(self.options.tspan);
                    end

                otherwise
                    error('Unrecognized solver');
            end

            if downsample
                self.nsys_output.a_hip.resample(downsample);
                self.nsys_output.a_isp.resample(downsample);
            end

            if any(burn)
                self.nsys_output.a_hip.burn(burn(1));
                self.nsys_output.a_isp.burn(burn(1));
                self.nsys_output.a_hip.rem_after(self.options.tspan);
                self.nsys_output.a_isp.rem_after(self.options.tspan);
            end
        end

        function extend(self,varargin)
            % Extend an existing run, discarding original results for memory efficiency
            % The ICs will be copied from the result, and the HIP and ISP quantities also copied
            % Unlike the base implementation, because self.param.init is changed, the runs are
            % reproducible. That is,
            %
            %   wc.run()
            %   wc.extend()
            %   wc.run()
            %
            % The last two runs should produce the same output (because the init variable was set by extend())
            % This won't work if the simulation sampling rate has been changed. If that is the case, wc.result()
            % should be resampled to the new sampling rate beforehand. 

            % Get the history length
            original_seed = self.param.noise_seed;
            original_init = self.param.init;
            original_units = self.units;

            try
                % Replace the initial conditions
                self.param.init = []; % Clear any existing history
                m = self.get_mex_inputs;
                history_length = length(m.problem.history.time);
                self.param.init = ra.TimeSeries(self.result.time(end-history_length-5:end),self.result.vals(end-history_length-5:end,:));

                % Copy plasticity variables
                for j = 1:self.n_units
                    self.units(j).cie = self.nsys_output.a_isp.vals(end,j);
                    self.units(j).E.S.mu = self.nsys_output.a_hip.vals(end,j);
                end

                self.param.noise_seed = self.param.noise_seed + 1; % Increment the seed

                self.run(varargin{:});

            catch ME % Restore original parameters in case of a failure
                self.param.init = original_init;
                self.units = original_units; 
                rethrow(ME)
            end
            
        end

        function mx_input = get_mex_inputs(self,tspan)
            if nargin < 2 || isempty(tspan) 
                tspan = self.options.tspan;
            end
                
            % Compute ivp, prop, sys, and other
            mx_input = struct();
            mx_input.problem = struct();
            mx_input.properties = struct();
            mx_input.system = struct();
            mx_input.other = struct();
            
            % SET PROP
            mx_input.properties.step.dt = self.options.tstep;
            mx_input.properties.error.abs_err = 1e-6;
            mx_input.properties.error.rel_err = 1e-3;
            
            % SET SYS
            mx_input.system.nodes = self.build_nodes(tspan);
            mx_input.system.edges = self.build_edges();
            mx_input.system.options = self.options;
            mx_input.system.options.enable_isp = self.param.enable_isp; % TODO: This should be inferred by nsys by checking for nonzero eta/rho

            if self.param.round_delays
                round_tstep = mx_input.properties.step.dt;
                for j = 1:length(mx_input.system.edges)
                    mx_input.system.edges(j).delay = round(mx_input.system.edges(j).delay/round_tstep)*round_tstep;
                end
            end

            % How many timesteps do the delays correspond to?
            d = [mx_input.system.edges.delay];
            if min(floor(d(d>0)/mx_input.properties.step.dt)) < 2 % With interpolation of delays, this is not as significant an issue, as long as it isn't smaller than 1 (otherwise would need to interpolate the solution)
                warning('Velocity is high enough that the delay is less than 2 timesteps. Decrease timestep or decrease velocity')
            end

            % SET IVP
            mx_input.problem.t_start = 0;
            mx_input.problem.t_end = tspan;
            
            % Set history size
            extend_to = (ceil(max(d)./mx_input.properties.step.dt)+3)*mx_input.properties.step.dt; % Include at least 3 in the history
            mx_input.problem.history.time = [-extend_to:mx_input.properties.step.dt:0].';
            mx_input.problem.history.vals = ones(length(mx_input.problem.history.time),self.n_states);

            if ~isempty(self.param.init) % User specified initial conditions
                if isnumeric(self.param.init)
                    if ~isvector(self.param.init) || length(self.param.init)~=self.n_states
                        error('Numeric initial conditions must be a vector of length n_states')
                    end
                    mx_input.problem.x_start = self.param.init(:);
                    mx_input.problem.history.vals = bsxfun(@times,mx_input.problem.history.vals,self.param.init(:).');

                elseif isa(self.param.init,'ra.TimeSeries')
                    % Interpret the last portion as the initial conditions

                    if ~(self.param.init.is_arithmetic)
                        error('Old result must be uniformly spaced');
                    end

                    if abs(mean(diff(self.param.init.time))-mx_input.properties.step.dt)>1e-7
                        keyboard
                        fprintf('History timestep: %.4e, Simulation timestep: %.4e\n',mean(diff(self.param.init.time)),mx_input.properties.step.dt);
                        error('Initial conditions must have the same time spacing as the new simulation');
                    end

                    if size(self.param.init.vals,1) < size(mx_input.problem.history.vals,1)
                        error('Old result is shorter than the required history size');
                    end

                    mx_input.problem.history.vals = self.param.init.vals(end-(size(mx_input.problem.history.vals,1)-1):end,:);

                else
                    error('Initial conditions (wc.param.init) should be a vector or a TimeSeries');
                end

            else % Default initial conditions 
                s = RandStream('mt19937ar');
                s.reset(self.param.init_seed);
                mx_input.problem.x_start = rand(s,1,self.n_states);
                mx_input.problem.history.vals = bsxfun(@times,mx_input.problem.history.vals,mx_input.problem.x_start );
                e_noise = arrayfun(@(x) x.E.P.sigma,self.units);
                i_noise = 2*arrayfun(@(x) x.I.P.sigma,self.units);
                fullnoise = bsxfun(@times,reshape([e_noise(:),i_noise(:)]',1,[]),randn(size(mx_input.problem.history.vals)));
                mx_input.problem.history.vals = mx_input.problem.history.vals + fullnoise;
                mx_input.problem.history.vals(mx_input.problem.history.vals<0) = 0; 
            end
           
        end

        function nodes = build_nodes(self,tspan)
            % Convert units into nsys nodes
            nodes(self.n_nodes) = self.units(end).I; % Put the last one in place to preallocate

            s = RandStream('mt19937ar');
            s.reset(self.param.noise_seed);
            t = 0:self.options.tstep:tspan;
            for j = 1:self.n_units
                nodes((j*2-1)) = self.units(j).E;
                nodes(j*2) = self.units(j).I;
                nodes((j*2-1)).P = ra.utils.input2timecourse(nodes((j*2-1)).P,t,s);
                nodes(j*2).P = ra.utils.input2timecourse(nodes(j*2).P,t,s);
            end
        end

        function edges = build_edges(self)
            n_edges = 4*self.n_units + self.network.n_edges;

            ptr = 1; % Index of this edge
            edges(n_edges) = struct('src',NaN,'dst',NaN,'delay',0,'coupling',0);

            function ptr = set_edge(ptr,src,dst,delay,coupling)
                if coupling ~= 0 % Disable this in future if we want connections that can start at 0 but increase in strength
                    edges(ptr) = struct('src',src,'dst',dst,'delay',delay,'coupling',coupling);
                    ptr = ptr + 1;
                end
            end

            % First, set local connections and delays
            for j = 1:self.n_units
                u = (j-1)*2; % E = u+1, I = u+2
                ptr = set_edge(ptr,u+1,u+1,0,self.units(j).cee);
                ptr = set_edge(ptr,u+1,u+2,0,self.units(j).cei);
                ptr = set_edge(ptr,u+2,u+1,0,self.units(j).cie);
                ptr = set_edge(ptr,u+2,u+2,0,self.units(j).cii);
            end
            
            % Next, set network connections and delays
            edgedata = table2array(self.network.edges);
            for j = 1:self.network.n_edges
                u = (edgedata(j,1)-1)*2; % REGION U, E1 = u+1, I1 = u+2
                v = (edgedata(j,2)-1)*2; % REGION V, E2 = v+1, I2 = v+2
                ptr = set_edge(ptr,u+1,v+1,edgedata(j,4)./self.param.velocity,self.param.coupling(1)*edgedata(j,3)); % E1 -> E2
                ptr = set_edge(ptr,u+1,v+2,edgedata(j,4)./self.param.velocity,self.param.coupling(2)*edgedata(j,3)); % E1 -> I2 
            end

            edges = edges(1:ptr-1); % Trim to the actual number of non-zero strength connections

        end
    end
end
