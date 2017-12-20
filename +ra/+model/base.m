classdef base < matlab.mixin.Copyable
    % Abstract base class for models that interface with NSYS
    % Models such as Wilson-Cowan or Robinson inherit from this

    % All models must have these, but an empty value is not sensible, so this is abstract
    properties (Abstract) 
        nodes_per_unit; % Number of populations per unit/brain region
        states_per_node; % Number of state variables for each population
        mx_helper; % Name of the mex wrapper to use when calling NSYS
    end

    % All models must have these, and empty values are sensible, so doesn't need to be abstract
    properties 
        param; % Parameters that characterize the system being simulated
        options; % Parameters that affect the numerics of the integration
        units; % Local parameters for individual units
        network; % Network object specifying regions and long range connnectivity
        version; % Version number for automatic upgrading in the future
        result; % Timeseries of the state variables (available after simulating)
        nsys_output; % Any other additional outputs from NSYS 
    end
       
    properties (Dependent)
        n_units; % Total number of units being simulated (from network graph)
        n_nodes; % Total number of nodes being simulated (n_units * nodes_per_unit)
        n_states; % Total number of state variables (n_nodes * states_per_node)
    end
    
    methods (Abstract)
        mx_input = get_mex_inputs(tspan); % Function to return mex inputs for nsys. tspan is optional (see run method below)
    end

    methods

        function extend(self,varargin)
            % Extend an existing run, discarding original results for memory efficiency
            % Suppose m = self.get_mex_inputs
            % If:
            %   - self.results is not empty
            %   - the results timestep is the same as m.properties.step
            %   - isfield(m.problem,'history') is true
            % Then:
            %   - The last part of self.results is used to replace m.history
            %   - The existing results will be discarded
            %   - The new results will be obtained, with new time bounds
            %
            % In practice, if you perform a run using a HistoryProblem, you can
            % change the model's tspan and then run extend() to perform a new run
            % with ICs obtained from the previous run, as long as you haven't changed
            % the timestep. All usual arguments to run() are supported
            %
            % This is a low-level implementation that does not provide run reproducibility
            % That is
            %
            %     base.run()
            %     base.extend()            
            %     base.run()
            %
            % The last two runs will be different, because they will have different initial conditions.
            % Classes may handle their ICs differently, so extending based on ICs should be implemented
            % by derived classes c.f. ra.model.WilsonCowan

            opt        = ra.utils.kwArgs(varargin{:});
            burn       = opt.get('burn',[0 0]); % Two elements corresponding to seconds to burn from the start and end
            downsample = opt.get('downsample',false); % downsample to a given sampling rate 
            m          = opt.get('mexcfg',self.get_mex_inputs(self.options.tspan + sum(burn) ));

            if isempty(self.result)
                error('Old result must be present to extend a run');
            end

            if ~(self.result.is_arithmetic)
                error('Old result must be uniformly spaced');
            end

            if ~(self.result.fs == 1/m.properties.step.dt)
                error('Old result must be at the same time spacing as the current step size');
            end

            if ~isfield(m.problem,'history')
                error('Run extending only works with HistoryProblem models at the moment');
            end

            if size(self.result.vals,1) < size(m.problem.history.vals,1)
                error('Old result is shorter than the required history size');
            end

            assert(size(self.result.vals,2)==size(m.problem.history.vals,2)); % Same number of state variables, unlikely to be different unless user has changed the network in the meantime...

            m.problem.x_start = self.result.vals(end,:);
            m.problem.history.vals = self.result.vals(end-size(m.problem.history.vals,1)+1:end,:); 

            self.run('burn',burn,'downsample',downsample,'mexcfg',m);
        
        end

        function run(self,varargin)
            % Run the simulation and populate the results
            opt        = ra.utils.kwArgs(varargin{:});
            burn       = opt.get('burn',[0 0]); % Two elements corresponding to seconds to burn from the start and end
            downsample = opt.get('downsample',false); % downsample to a given sampling rate 
            m          = opt.get('mexcfg',self.get_mex_inputs(self.options.tspan + sum(burn) ));

            % If self.get_mex_inputs returned a structure containing a HistoryProblem, make sure it is valid
            if isfield(m.system.options,'fixed_step_size') && m.system.options.fixed_step_size && isfield(m.problem,'history') && max(abs((diff(m.problem.history.time)./m.system.options.tstep)-1)) > 1e-7
                error('If using HistoryProblem and fixed_step_size, the history timestep must be the same as the integration timestep');
            end

            self.nsys_output = mx.(self.mx_helper)(m);
            self.result = ra.TimeSeries(self.nsys_output.solution.time,self.nsys_output.solution.vals);
            self.nsys_output = rmfield(self.nsys_output,'solution'); % Save space by discarding duplicate information
            self.result.rem_before(0);
            
            if downsample
                self.result.resample(downsample);
            end

            if any(burn)
                self.result.burn(burn(1));
                self.result.rem_after(self.options.tspan);
            end

        end

        function n = get.n_units(self)
            % Total number of units being simulated (from network graph)
            n = length(self.units); 
        end

        function n = get.n_nodes(self)
            % Total number of nodes being simulated
            n = self.n_units*self.nodes_per_unit;
        end

        function n = get.n_states(self)
            % Total number of state variables
            n = self.n_nodes*self.states_per_node;
        end

        function clear(self)
            % Delete results and nsys_output to reduce the size of the object
            self.result = [];
            self.nsys_output = [];
        end

        function plot(self)
            % Basic plot functionality
            plot(self.result);
        end

        function analyse(self)
            % Call a basic standard analysis
            self.plot()
        end
        
    end

    methods(Access = protected)
        function new = copyElement(self)
            new = copyElement@matlab.mixin.Copyable(self); % Shallow copy things by default
            % Deep copy results
            if ~isempty(self.result)
                new.result = self.result.copy();
            end
        end
    end
end
