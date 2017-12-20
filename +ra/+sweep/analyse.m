classdef analyse < Pepper
	% This class internally stores the unit, network, velocity, and sweeps over them
	% Intended usage is to make a file like 'sweep_deco_2.m' containing
	%
	% function p = sweep_deco_2()
	%	  u = ra.unit.deco;
	%	  u.E.P.mu = 0.35;
	%     p = ra.sweep.wc_pepper(ra.network.line(2),u,linspace(1,80,3),linspace(0,0.5,3))
	%
	% And then doing './submit.py sweep_deco_2.m'
	
	properties
		fname
	end

	methods(Static)
		function out = postprocess(inputs,results)
			% In the postprocess function, extract from cell and reshape
			out_temp = postprocess@Pepper(inputs,results);
			f = fields(out_temp);
			for j = 1:length(f)
				if ischar(out_temp(1).(f{j}))
					out.(f{j}) = reshape({out_temp.(f{j})},size(inputs.simdata.inputs.c,1),size(inputs.simdata.inputs.c,2),size(inputs.simdata.inputs.c,3));
				else
					try
	                    out.(f{j}) = reshape([out_temp.(f{j})],size(inputs.simdata.inputs.c,1),size(inputs.simdata.inputs.c,2),size(inputs.simdata.inputs.c,3));
	                catch
	                    out.(f{j}) = reshape({out_temp.(f{j})},size(inputs.simdata.inputs.c,1),size(inputs.simdata.inputs.c,2),size(inputs.simdata.inputs.c,3));
	                end
	            end
			end
		end
	end

	methods 

		function self = analyse(fname)
			% analyse takes in the filename of results.mat - or it could be results_inputs.mat
			self.fname = fname
		end

		function n_workers = assemble(self)
			% Move the analysis back into its own file
			n_workers = assemble@Pepper(self);
			movefile('results.mat',strrep(self.fname,'results.mat','analysis.mat'));
		end

		function [inputs,task_ids] = get_inputs(self)
			inputs.simdata = load(self.fname);
			%task_ids = 1:numel(inputs.simdata.inputs.task_ids); % Task ids specify which jobs are to be run, matches the original number of tasks
			task_ids = 1:numel(inputs.simdata.inputs.c);
		end

		function out = do_work(self,task_id,inputs)
			original_inputs = inputs.simdata.inputs;

			d = load(fullfile(fileparts(self.fname),sprintf('results_part_%d.mat',task_id))); % Load results_part*.mat starting with path to results.mat
			wc = d.out.wc;
			wc.result.cast('double');
			isp_ts = d.out.isp_ts;
			out.fname = sprintf('results_part_%d.mat',task_id);
			
			% If we have done a long simulation run from simulate_plasticity.m, then the output would not
			% have been trimmed as part of the integration e.g. wc.run('burn',[5 5]). Thus we need to make sure
			% it is done here
			if wc.result.tspan > 500
				wc.result.trim([15 10]); % Discard from the start and end to account for the downsampling edge effects
			end

			excitatory = wc.excitatory;
			inhibitory = wc.inhibitory;

			% Fit cie and node strength as a measure of EI balance
			str = wc.network.strength.';
			p = polyfit(str,wc.cie,1)
			cv = polyval(p,str);
			yresid = wc.cie - cv;
			SSresid = sum(yresid.^2);
			SStotal = (length(cv)-1) * var(wc.cie);
			out.ei_rsquared = 1 - SSresid/SStotal;
			out.ei_correlation = corr(str,wc.cie);
			
			% Compute spectral variability - no overlap, if they overlap we see chirps where there were none
		    window_length = round(1*excitatory.fs);
		    overlap = round(0*window_length);
		    bw = nan(excitatory.n_signals,1);
		    log_bw = nan(excitatory.n_signals,1);
		    for j = excitatory.n_signals:-1:1
		    	[~,f,time,P(:,:,j)] = spectrogram( detrend(excitatory.vals(:,j)), window_length,overlap, window_length, excitatory.fs );
		    end
			fvar = std(P,[],2); % Standard deviation over time
			fmean = mean(fvar,3); % Mean variability over channels, for each frequency
			out.power_std = trapz(f,fmean); % Integrate over frequency

			% Compute stats from single channel spectra
			% Note - peak frequency must be > 1 (to prevent very low freq. components from dominating)
		    bw = nan(excitatory.n_signals,1);
		    log_bw = nan(excitatory.n_signals,1);
		    for j = excitatory.n_signals:-1:1
		    	% Compute bandwidth using a longer window to get better spectral resolution
		    	[output_P(:,j),out.spec_f] = pwelch(detrend(excitatory.vals(:,1)),round(5*excitatory.fs),[],[],excitatory.fs); % Get PSD
	    		channel_f = out.spec_f;
	    		channel_P = output_P(:,j);
	    		bw(j) = mean(channel_P>(0.1*max(channel_P))); % What fraction of the spectrum has at least 10% of the power
	    		log_bw(j) = mean(log(channel_P)>(0.1*max(log(channel_P)))); % What fraction of the spectrum has at least 10% of the power
		   	
	    		channel_P = channel_P(channel_f>1);
	    		channel_f = channel_f(channel_f>1);
	    		[~,idx] = max(channel_P);
	    		maxp_gt1(j) = channel_f(idx);
		    end
		    out.power_bw = mean(bw);
		    out.power_bw_log = mean(log_bw);
		    out.f_max = mean(maxp_gt1);
		    out.f_std = std(maxp_gt1);
		    out.spec_P = mean(output_P,2); % Save the power spectrum to make plotting quicker

			% Do HMM fits
			% [~,out.hmm_stats] = ra.hmm_regression(wc);
			% [~,out.hmm_orth_stats] = ra.hmm_regression(wc,[],true);

			% Compute the oscillatory amplitude in each region
			osc_amp = std(excitatory.vals);
			out.amp_ratio = max(osc_amp)./min(osc_amp);

			[~,~,spike_rate] = ra.spike_counter(wc);
			out.spike_rate = mean(spike_rate);

			% Overall mean network activity 
			out.excitatory_mean = mean(excitatory.vals(:));
			out.inhibitory_mean = mean(inhibitory.vals(:));

			% Compute IWE activity
			iwe = wc.iwe;
			out.iwe_discrepancy = max(iwe)-min(iwe);
			out.iwe_plastic_discrepancy = max(abs(iwe-wc.units(1).I.isp_target));

			% Std of time-averaged activity across ROIs
			out.excitatory_std = std(mean(excitatory.vals));
			out.inhibitory_std = std(mean(inhibitory.vals));

			% Discrepancy in time-averaged activity across ROIs 
			out.excitatory_discrepancy = max(mean(excitatory.vals))-min(mean(excitatory.vals));
			out.inhibitory_discrepancy = max(mean(inhibitory.vals))-min(mean(inhibitory.vals));

			% Mean oscillatory amplitude
			out.excitatory_amp = mean(std(excitatory.vals));
			out.inhibitory_amp = mean(std(inhibitory.vals));

			% Variability in oscillatory amplitude
			out.excitatory_amp_std = std(std(excitatory.vals));
			out.inhibitory_amp_std = std(std(inhibitory.vals));

			% Discrepancy in oscillatory amplitude
			out.excitatory_amp_discrepancy = max(std(excitatory.vals))-min(std(excitatory.vals));
			out.inhibitory_amp_discrepancy = max(std(inhibitory.vals))-min(std(inhibitory.vals));

			% Overall cei
			out.cie_mean = mean([wc.units.cie]);
			out.cie_std = std([wc.units.cie]);

			% Load in the dfilt and non-dfilt data
			data = load(fullfile(startup.get_rootdir,'mark_data','mark_final_connectivity_alpha.mat')); % Could use dfilt instead

			% Rescale simulation first
			ts = wc.excitatory.copy;
			out.model_freqs = data.band_frequency; % Use the same bands as the data
			% [~,idx] = min(abs(model_freqs-10)); % Find which model frequency index corresponds to the alpha band
			% [ts,out.f_adjust_factor] = ra.scale_frequency(ts,10,[0.5 2]); % Rescale frequency
			out.f_adjust_factor = 1;
			[~,Ph] = ts.envelope('filter',data.band_frequency); % Don't downsample
			order_ts = abs(mean(exp(1i*Ph.vals),2));
			out.synchrony = mean(order_ts);
			out.metastability = std(order_ts);
					
			% We are only looking at alpha band now, so just set it directly
			out.alpha_synchrony = out.synchrony;
			out.alpha_metastability = out.metastability;

			[out.max_synchrony] = NaN;
			out.max_synchrony_f = NaN;

			[out.max_metastability] = NaN;
			out.max_metastability_f = NaN;

			% Goodness of fit, WITHOUT rescaling the peak frequency
			[tmp,~,out.alpha_conn_aec] = ra.analysis.cmat_comparison_reduced(wc,'aec','orthogonalize',true,'data',data);
			out.alpha_conn_aec_matrix = tmp;
			[tmp,~,out.alpha_conn_plv] = ra.analysis.cmat_comparison_reduced(wc,'plv','orthogonalize',true,'data',data);
			out.alpha_conn_plv_matrix = tmp;
			[tmp,~,out.alpha_conn_pli] = ra.analysis.cmat_comparison_reduced(wc,'pli','orthogonalize',true,'data',data);
			out.alpha_conn_pli_matrix = tmp;
			[tmp,~,out.alpha_conn_no_orth_aec] = ra.analysis.cmat_comparison_reduced(wc,'aec','orthogonalize',false,'data',data);
			out.alpha_conn_no_orth_aec_matrix = tmp;
			[tmp,~,out.alpha_conn_no_orth_plv] = ra.analysis.cmat_comparison_reduced(wc,'plv','orthogonalize',false,'data',data);
			out.alpha_conn_no_orth_plv_matrix = tmp;
			[tmp,~,out.alpha_conn_no_orth_pli] = ra.analysis.cmat_comparison_reduced(wc,'pli','orthogonalize',false,'data',data);
			out.alpha_conn_no_orth_pli_matrix = tmp;
			
			out.raw_correlation = corr(ts.vals);
			tsf = ts.filter(data.band_frequency);
			out.raw_correlation_alpha = corr(tsf.vals);

			% Need this return below to work with the non-ISP runs
			if isempty(isp_ts) || isp_ts.tspan < 200
				return
			end

			isp_ts.burn(1000); % Skip the fast plasticity parts
			out.max_isp_variability = max(std(isp_ts));
			out.mean_isp_variability = mean(std(isp_ts));

		end

	end

end
