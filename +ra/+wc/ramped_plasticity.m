function [wc1,isp_ts,activity_ts] = ramped_plasticity(wc,tspan,n_runs,plasticity_rates)
	% Take in a WC, clear it, then run plasticity
	% Return an unsimulated WC ready to run with wc.run

	% Do long plasticity run

	% plasticity_rates first column stores time from which the rate takes effect, second column stores the rate itself
	if nargin < 4 || isempty(plasticity_rates) 
		plasticity_rates = [0 0.4; 500 0.1; 1000 0.05];
	end
	
	assert(all(mod(plasticity_rates(:,1)/tspan,1)==0),'Plasticity changeover times must be integer multiples of the stride time');

	t_count = tic;

	fprintf('Plasticity starting initial run, t=0\n');

	wc.options.tspan = tspan; % Go in increments of this 
	wc.clear();

	set_rate(wc,plasticity_rates,0);
	wc.run();

	fprintf('Finished initial run, t=%.2f\n',toc(t_count));

	% Set up low frequency output records
	low_fs = 1; % Downsample the activity over the full timeseries onto this scale
	isp_track = struct;
	activity_track = struct;
	isp_track = ra.utils.downsample_track(isp_track,wc.isp,low_fs);
	activity_track = ra.utils.downsample_track(activity_track,wc.result,low_fs);

	fprintf('Finished initial downsample, t=%.2f\n',toc(t_count));

	for j = 2:n_runs % The first run has already been completed

		fprintf('Plasticity iteration %d (t=%.2f)\n',j,toc(t_count))

		set_rate(wc,plasticity_rates,(j-1)*tspan); % On the first iteration, the start time is 1*tspan
		wc.extend();
		isp_track = ra.utils.downsample_track(isp_track,wc.isp,low_fs);
		activity_track = ra.utils.downsample_track(activity_track,wc.result,low_fs);
	end

	isp_ts = isp_track.output;
	activity_ts = activity_track.output;

	% Turn plasticity off and do a 1s extend so that the output WC is ready
	% to simulate with ISP off and with the correct initial conditions
	fprintf('Commencing final run, t=%.2f\n',toc(t_count));

	wc1 = wc.copy();
	for j = 1:wc1.n_units
		wc1.units(j).I.isp_rate = 0;
		wc1.units(j).E.hip_rate = 0;
	    wc1.units(j).cie = wc1.nsys_output.a_isp.vals(end,j);
	    wc1.units(j).E.S.mu = wc1.nsys_output.a_hip.vals(end,j);
	end
	wc1.options.tspan = 1;
	wc1.extend(); % Downsampling to 300Hz to hopefully reduce file size a bit further. Also comparable to MEG 
	wc1.clear();

end

function set_rate(wc,plasticity_rates,t)
	% Take in wc, plasticity rates, and current time
	% Set the plasticity rate accordingly
	idx = find(t>=plasticity_rates(:,1),1,'last');
	for j = 1:length(wc.units)
		wc.units(j).I.isp_rate = plasticity_rates(idx,2); 
	end
end

