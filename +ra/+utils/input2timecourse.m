function tc = input2timecourse( in, time , rstream )
%
% Transform a stimulus specified as a structure to a timecourse (structure with fields 'time' and 'vals').
% JH (2016?)

    % if the input is already a time-course, don't process it
    if all(isfield( in, {'time','vals'} ))
        tc = in; return;
    end
    
    % otherwise, create a timecourse
    tc.time = time(:);
    dt = time(2)-time(1);
    nt = numel(time);

    % Use optional RNG stream if provided with one
    if nargin < 3 || isempty(rstream) 
        rstream = RandStream.getGlobalStream;
    end

    switch lower(in.name)
        
        case {'const','constant'}
            tc.time = [tc.time(1) (tc.time(1)+tc.time(end))/2 tc.time(end)];
            tc.vals = in.value * ones(3,1);
            
        case {'gaussian','gp'}
            if in.sigma == 0
                tc = ra.utils.input2timecourse(struct('name','const','value',in.mu),time);
            else    
                tc.vals = in.mu + in.sigma * randn(rstream,nt,1);
            end
            
        case {'gaussian_tapered'}

            if in.sigma == 0
                tc = network.helper.input2timecourse(struct('name','const','value',in.mu),time);
            else   
                taper_time = tc.time > in.taper(1) & tc.time < in.taper(2);
                taper = +(tc.time<=in.taper(1));
                taper(taper_time) = linspace(1,0,sum(taper_time)); 
                taper = taper.*in.sigma;

                if isfield(in,'sigma_final')
                    assert(in.sigma_final < in.sigma);
                    taper(taper<in.sigma_final) = in.sigma_final;
                end
                
                tc.vals = in.mu + taper.*randn(rstream,nt,1);
            end

        case {'wiener','rw'}
            tc.vals = [0;cumsum( dk.get_field(in,'sigma',1)*sqrt(dt)*randn(rstream,nt-1,1) )];
            
        case {'colored'}
            if in.sigma == 0
                tc = network.helper.input2timecourse(struct('name','const','value',in.mu),time);
            else    
                x = dsp.ColoredNoise('Color',in.color,'SamplesPerFrame',n,'NumChannels',1);
                tc.vals = in.mu + in.sigma * x.step();
            end
            
        case {'brw'}
            if in.sigma == 0
                tc = network.helper.input2timecourse(struct('name','const','value',in.mu),time);
            else
                tc.vals = in.mu + in.sigma * jh.algo.brw( nt, 1, dt, in.eta, in.alpha, dk.get_field(in,'gamma',0) );
            end

        case {'harmonic','sinusoid'}
            x = 2*pi*in.omega*tc.time + dk.get_field(in,'phi',0);
            tc.vals = in.avg + in.amp * sin(x);
            
        otherwise
            error('Unknown input type: "%s".',in.name);
        
    end

end