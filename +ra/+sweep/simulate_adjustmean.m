classdef simulate_adjustmean < Pepper
	% The purpose of this file is to mix existing results with new runs
	% Historically this was used to go from deco_rk4_cortical to deco_rk4_hd
	
	properties
		unit % The unit being tested
		network % The network being tested
		velocity % An array of velocities
		coupling % An array of couplings
		cie_file
	end

	methods 

		function self = simulate_adjustmean(network,unit,velocity,coupling,cie_file)
			self.unit = unit;
			self.network = network;
			self.velocity = velocity;
			self.coupling = coupling;
			self.cie_file = cie_file;
		end

		function n_workers = assemble(self)
			% Save only the inputs and the number of workers
			% Since the results are split across
			% results.mat and result_part_**.mat
			[inputs,task_ids] = self.get_inputs();
			d = load('worker_files/worker_output_1.mat');
			n_workers = d.n_workers;
			save('results.mat','inputs','n_workers','task_ids');

			% Automatically run the analysis on the first run, with the same number of workers as the original
			[rootdir,run_name] = fileparts(pwd());
			cd('..');
			fprintf('SUBMITTING: source ~/.bash_profile && ssh jalapeno "cd %s && submit %d analyser.m --argstring \\"(''%s'')\\""\n',rootdir,n_workers,run_name);
			system(sprintf('ssh jalapeno "source ~/.bash_profile && cd %s && submit %d analyser.m --argstring \\"(''%s'')\\""',rootdir,n_workers,run_name));
		end

		function [inputs,task_ids] = get_inputs(self)
			[inputs.c,inputs.v] = meshgrid(self.coupling,self.velocity);
			inputs.velocity = self.velocity;
			inputs.coupling = self.coupling;
			task_ids = reshape(1:numel(inputs.v),size(inputs.v));			
		end

		function out = do_work(self,task_id,inputs)
			% Check if this case has already been run
			out = struct();

			this_run_fname = sprintf('results_part_%d.mat',task_id);

			if exist(this_run_fname)
				try
					load(this_run_fname)
					fprintf('Successfully loaded %d - skipping\n',task_id)
					return
				end
			end

			cie_vals = load(self.cie_file);
			u = self.unit;
			u.I.isp_rate = 0; % Make sure ISP rate is zero
			u.cie = cie_vals.cie_mean(task_id);
			cie_vals.cie_mean(task_id)
			
			% Initialize the model and define integration settings
			wc = ra.model.WilsonCowan(u,self.network);
			wc.param.coupling = [inputs.c(task_id) 0];
			wc.param.velocity = inputs.v(task_id);
			wc.options.tstep = 1e-4; % Just barely OK for a 1ms mean delay. Should be fine for 2ms or more. 
			wc.mx_helper = 'wc_network_rk4';
			wc.param.round_delays = false;

			% Do the output run
			wc = ra.wc.long_run(wc,25,21,300); % PROPER RUN
			% wc = ra.wc.long_run(wc,25,3,300); % FAST RUN

			save(this_run_fname,'wc','-v7.3'); 

			wc.result.vals = single(wc.result.vals); % Save as a single
			wc.nsys_output = []; % No need to retain HIP and ISP variables

			out.wc = wc;
			out.isp_ts = [];
			out.activity_ts = [];

			save(this_run_fname,'out','-v7.3');

			out = struct('completed',true); % Skip saving outputs
		end

	end

end

