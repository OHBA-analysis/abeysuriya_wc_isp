classdef TimeSeries < matlab.mixin.Copyable
    % Class to store timeseries data consisting of time, and one or more variables
    % Time is assumed to be in seconds
    % Rows of values correspond to observations
    % JH and RA
    
    properties
        time
        vals
    end

    properties (Dependent)
        n_times
        n_signals
        tspan
    end
    
    properties (Dependent,Hidden)
        nt
        ns
    end

    methods
        
        function self = TimeSeries(t,v)
            % Time can be a sampling rate (scalar) or time values (vector)
            % V should be a matrix where the number of rows or columns matches the number of time points
            % If time is a scalar, then rows of V are assumed to be time points
            if nargin == 0
                return
            end

            if isa(t,'meeg')
                v = t(:,:,:).';
                t = t.time;
            end

            if any(~isreal(t) | ~isfinite(t))
                error('Time values must be finite and real')
            end

            if ~isvector(t)
                error('Time must be a vector');
            end

            % Interpret as a sampling rate - assume that v is oriented with rows corresponding to timepoints
            if length(t) == 1 
                t = 0:1/t:1/t*(size(v,1)-1);
            end

            self.time = t(:);

            assert( isnumeric(v) && ismatrix(v) && any(size(v) == size(self.time,1)), 'Values should be n_times x n_signals.' );

            if ndims(v) > 2
                error('V must be a 2D matrix');
            end

            if size(v,1) ~= length(t)
                v = v.';
            end
            
            self.time = t(:);
            self.vals = v;
        end
                
        % Dimensions
        function n = get.nt(self), n = size(self.time,1); end % num of timesteps
        function n = get.ns(self), n = size(self.vals,2); end % num of signals
        function n = get.n_times  (self), n = self.nt; end
        function n = get.n_signals(self), n = self.ns; end
        
        function Dt = get.tspan(self)
            Dt = self.time(end) - self.time(1);
        end
                
        % Arithmetic if all timesteps are almost equal
        function yes = is_arithmetic(self,thresh)
            if nargin < 2, thresh = 5e-5; end
            yes = std(diff( self.time )) <= thresh;
        end
        
        % Interpolate at given timepoints
        function ts = interpolate(self,query_t,method)
            if nargin < 3, method = 'pchip'; end
            
            query_t  = query_t(:);
            new_vals = interp1( self.time, self.vals, query_t, method );
            if nargout == 0
                self.vals = new_vals;
                self.time = query_t;
            else
                ts = ra.TimeSeries( query_t, new_vals );
            end
        end
        
        % Resample at a given sampling frequency
        function ts = resample(self,fs,method)
            if nargin < 3, method = 'pchip'; end

            meanvals = nanmean(self.vals);

            if self.n_times > 1e5 % For large arrays, split the task to save memory
                [new_vals(:,self.n_signals),new_time] = resample(self.vals(:,self.n_signals)-meanvals(self.n_signals), self.time, fs, method);
                for j = 1:self.n_signals-1
                    [new_vals(:,j),new_time] = resample(self.vals(:,j)-meanvals(j), self.time, fs, method);
                end
            else
                [new_vals,new_time] = resample( bsxfun(@minus,self.vals,meanvals), self.time, fs, method );
            end

            new_vals = bsxfun(@plus,new_vals,meanvals);

            if nargout == 0
                self.vals = new_vals;
                self.time = new_time;
            else
                ts = ra.TimeSeries( new_time, new_vals );
            end
        end
        
        function n = fs(self)
            if self.is_arithmetic
                n = 1/(self.time(2)-self.time(1));
            else
                n = NaN;
            end
        end
                
        % Mean and sdev
        function m = mean(self,varargin), m = mean( self.vals, varargin{:} ); end
        function s = std(self,varargin), s = std ( self.vals, varargin{:} ); end
        
        % Time masks
        function m = tmask_lt(self,cut), m = round(self.time,8) < cut; end
        function m = tmask_gt(self,cut), m = round(self.time,8) > cut; end
        
        function m = tmask_leq(self,cut), m = round(self.time,8) <= cut; end
        function m = tmask_geq(self,cut), m = round(self.time,8) >= cut; end
        
        function ts = burn(self,tval)
            if ~self.is_arithmetic
                error('Not checked properly if input TS is not arithmetic');
            end

            % If arithmetic, all times should be integer multiples of the step size
            dt = self.time(2)-self.time(1);
            self.time = dt*round(self.time./dt);

            if nargout == 0
                self.rem_before(tval);
                self.time = self.time - tval;
                self.time = self.time-self.time(1);
            else
                ts = self.rem_before(tval);
                ts.time = ts.time - tval;
                ts.time = ts.time - ts.time(1);
            end
        end

        function ts = remove_times(self,sel)
            if islogical(sel)
                sel = ~sel;
            else
                sel = setdiff( 1:self.nt, sel );
            end
            if nargout == 0
                self.select_times(sel);
            else
                ts = self.select_times(sel);
            end
        end
      
        % Remove timepoints after/before a certain time
        function ts = rem_before(self,tval)
            m = self.tmask_lt(tval);
            if nargout == 0
                self.remove_times(m);
            else
                ts = self.remove_times(m);
            end
        end
        function ts = rem_after(self,tval)
            m = self.tmask_gt(tval);
            if nargout == 0
                self.remove_times(m);
            else
                ts = self.remove_times(m);
            end
        end
        
        % Select a subset of signals or timepoints (in-place if no output)
        function ts = select_signals(self,sel)
            if nargout == 0
                self.vals = self.vals(:,sel);
            else
                ts = ra.TimeSeries( self.time, self.vals(:,sel) );
            end
        end
        
        function ts = select_times(self,sel)
            if nargout == 0
                self.time = self.time(sel);
                self.vals = self.vals(sel,:);
            else
                ts = ra.TimeSeries( self.time(sel), self.vals(sel,:) );
            end
        end
        
        % Remove a subset of signals or timepoints (in-place if no output)
        function ts = remove_signals(self,sel)
            if islogical(sel)
                sel = ~sel;
            else
                sel = setdiff( 1:self.ns, sel );
            end

            if nargout == 0
                self.select_signals(sel);
            else
                ts = self.select_signals(sel);
            end
        end

        function cast(self,type)
            self.time = cast( self.time, type );
            self.vals = cast( self.vals, type );
        end

    end
    
    
    %-------------------------
    % Value-changing methods
    %-------------------------
    methods
                
        function ts = rereference(self,ref)
            % Re-reference the cross-signal average with respect to a subset of signals.
            if nargin < 2
                ref = 1:self.ns;
            end

            ref = ref(:).';
            mu  = mean( self.vals(:,ref), 2 );
            v   = bsxfun( @minus, self.vals, mu );

            if nargout == 0, 
                self.vals = v;
            else            
                ts = ra.TimeSeries( self.time, v );
            end
        end
        
        function ts = orthogonalize(self,method)
            % Orthogonalize using ROInets.remove_source_leakage
            if nargin < 2 || isempty(method) 
                method = 'symmetric';
            end
                
            v = ROInets.remove_source_leakage(self.vals.', method).';

            if nargout == 0
                self.vals = v;
            else
                ts = ra.TimeSeries(self.time,v);
            end
        end

        function varargout = filter(self,fband)
            % Apply a filter using FT_PREPROC functions
            %
            % bandpass - [8 12]
            % lowpass - [0 5]
            % highpass - [10 inf]

            if ischar(fband)
                f = nsl.sig.priv.fbander_parse(fband);
                fband = f.freq;
            end
            dat = self.vals.';

            % Filter order presently corresponds to fieldtrip defaults
            if fband(1) == 0 
                dat2 = ft_preproc_lowpassfilter(dat, self.fs, fband(2), 6, 'but','twopass','no');
            elseif isinf(fband(2))
                dat2 = ft_preproc_highpassfilter(dat, self.fs, fband(1), 6, 'but','twopass','no');
            else
                dat2 = ft_preproc_bandpassfilter(dat, self.fs, fband, 4, 'but','twopass','no');
            end

            if nargout == 0
                self.vals = dat2.';
            else
                varargout{1} = ra.TimeSeries(self.time,dat2.');
            end
        end

        % Compute the numerical derivative of the signal
        function ts = derivative(self,h)
            dv = diff(self.vals);
            dt = diff(self.time);
            dvdt = dv./dt;
            t = self.time(1:end-1)+dt/2;

            if nargout == 0
                self.time = t;
                self.vals = dvdt;
            else
                ts = ra.TimeSeries(t,dvdt);
            end
        end

        
        function [Hen,Phase] = envelope(self,varargin)
            % Return Hilbert envelope and corresponding phases
            %
            % INPUTS
            % - filter : specify range of frequencies to filter at prior to enveloping (e.g. [8 12], default: no filter)
            % - rescale : Rescale envelopes to unit mean and unit variance (default: false)
            % - downsample : Resample the envelopes and phase timeseries at this frequency, also trim 1s on either end 
            %                (default: no downsampling). Set scalar to use same frequency for both. Set NaN to skip
            % - orthogonalize: Orthogonalize time series after filtering but before enveloping (default: false)

            arg = inputParser;
            arg.addParameter('filter',[]);
            arg.addParameter('rescale',false);
            arg.addParameter('downsample',[]);
            arg.addParameter('orthogonalize',false);
            arg.parse(varargin{:});

            ts = self.copy();

            if ~isempty(arg.Results.filter)
                ts.filter(arg.Results.filter);
            end

            if arg.Results.orthogonalize
                ts.orthogonalize;
            end

            analytic = hilbert(bsxfun( @minus, ts.vals, ts.mean ));
            Hen = ra.TimeSeries(ts.time,  abs(analytic));
            Phase = ra.TimeSeries(ts.time,  unwrap(angle(analytic)));

            if ~isempty(arg.Results.downsample)
                ds = arg.Results.downsample;
                if length(ds) == 1
                    ds(2) = ds(1);
                end

                if isfinite(ds(1))
                    Hen = Hen.resample(ds(1));
                    Hen.burn(1);
                    Hen.rem_after(Hen.tspan-1);
                end

                if isfinite(ds(2))
                    Phase = Phase.resample(ds(2));
                    Phase.burn(1);
                    Phase.rem_after(Phase.tspan-1);
                end

            end

            if arg.Results.rescale
                Hen.vals = bsxfun(@minus,Hen.vals,mean(Hen.vals)); 
                Hen.vals = bsxfun(@rdivide,Hen.vals,std(Hen.vals)); 
            end

            if nargout == 0
                self.time = Hen.time;
                self.vals = Hen.vals;
            end

        end

        % General plot - timecourse
        function h = plot(self,varargin)
            h = plot( self.time', self.vals', varargin{:} );
            h = h(1);
            ax = get(h,'Parent');
            set(ax,'XLim',[self.time(1) self.time(end)]);
            xlabel(ax,'Time (s)')
            ylabel(ax,'Signal value')
        end
                
        % Image plot (for lots of signals)
        function tc_image(self,varargin)
            img = imagesc( self.time', 1:self.ns, self.vals', varargin{:} ); 
            ax = get(img,'Parent');
            xlabel(ax,'Time (s)');
            ylabel(ax,'Signal');
            cbar = colorbar(ax);
            ylabel(cbar,'Signal value');
            set(ax,'YDir','normal')
        end

        function [f,P] = tc_spectrum(self,idx)
            if nargin < 2 || isempty(idx) 
                idx = 1:self.n_signals;
            end
            
            legend_entries = {};
            figure
            for j = 1:length(idx)
                [P,f] = pwelch(self.select_signals(idx).vals-mean(self.select_signals(idx).vals),[],[],[],self.fs);
                loglog(f,P);
                hold on
                legend_entries{j} = sprintf('Signal %d',idx(j));
            end

            xlabel('Frequency (Hz)');
            ylabel('Power spectral density ([data units]^2 s^{-1})')
            legend(legend_entries)

        end

        function ts = trim(self,t)
            % Remove from start and end 
            %
            % Inputs
            % t - Scalar equivalent to [t t]. Remove t(1) from start and t(2) from end

            if length(t) == 1
                t(2) = t(1);
            end

            if nargout == 0
                self.burn(t(1));
                self.rem_after(self.tspan-t(2));
            else 
                ts = self.burn(t(1));
                ts.rem_after(ts.tspan-t(2));
            end
        end

        function tsw = window(self,winlen,overlap)
            % Window the timeseries into windows
            % INPUTS
            % - winlen : Window length in seconds (e.g. 5)
            % - overlap : Fraction of overlap when sliding windows (e.g. 0.5)
            % OUTPUTS
            % - Array of time series objects corresponding to each of the windows

            winlen_s = round(winlen*self.fs);
            stride = winlen_s - round(winlen_s*overlap);
            assert(winlen_s > 0,'Window length is too small compared to sampling rate');
            start_idx = 1:stride:self.n_times;
            stop_idx = start_idx + winlen_s;
            flt = stop_idx <= self.n_times;
            start_idx = start_idx(flt);
            stop_idx = stop_idx(flt);

            for j = 1:length(start_idx)
                tsw(j) = self.select_times(start_idx(j):stop_idx(j));
            end

        end

    end
    
end
