function wc = long_run(wc,tspan,n_runs,low_fs)
	% Do a long simulation and downsample result
	
	wc.options.tspan = tspan; % Go in chunks of this size
	wc.clear()
	activity_track = struct;

	t_start = tic;
	fprintf('Commencing initial run...\n')
	wc.run();
	fprintf('Initial run finished, t=%.2f\n',toc(t_start));

	activity_track = ra.utils.downsample_track(activity_track,wc.result,low_fs);

	for j = 2:n_runs
		fprintf('Starting run %d, t=%.2f\n',j,toc(t_start));
		wc.extend();
		fprintf('Finished run %d, t=%.2f\n',j,toc(t_start));
		activity_track = ra.utils.downsample_track(activity_track,wc.result,low_fs);
	end

	wc.result = activity_track.output;
