function d = downsample_track(d,new_ts,low_fs)
	% This function operates on a struct to downsample a timeseries that is incrementally extended
	% In the downsampling, both the start and end of the timeseries are affected by edge effects
	% However, if we consider
	% t1 = [0 10], t2 = [10 20], t3 = [20 30]
	% First, we cache t1 and t2. We stitch them into a single timeseries e.g.
	% p1 = [0 20]
	% Downsampling, we might have clean data from [1 19]
	% Next, t3 arrives, and we stitch to obtain
	% p2 = [10 30]
	% Downsampling, we have clean data from [11 29]
	% So we stitch the [1 19] and [11 29] together, because they match up in the clean parts e.g. 
	% t=15. So we take [0 15] from the first timeseries, and [15 30] from the second timeseries
	% Now the discontinuity in the middle has been removed

	% d is a struct with fields output (downsampled timeseries) and input
	if isempty(fields(d))
		d.prev = new_ts;
		d.output = d.prev.resample(low_fs);
		d.iterations = 1;
	else
		d.iterations = d.iterations + 1;

		i3 = append(d.prev,new_ts); 
		low_fs = d.output.fs;
		i3.resample(low_fs); % This spans i1 and i2

		% Move the new sample to start halfway through the old one
		% Note we have assumed that i3.time(1) == 0
		assert(i3.time(1) == 0) % If this is not zero, the indexing below won't be correct

		if d.iterations > 2 % If this is the second iteration, i3 really does start at zero. Otherwise, i3 starts d.prev.tspan beforehand
			cutover_time = d.output.time(end)-d.prev.tspan;
			i3.time = i3.time + cutover_time;
		end

		cutoff_time = d.output.tspan-5;
		new_time = [d.output.time(d.output.time <= cutoff_time);i3.time(i3.time > cutoff_time) ];
		new_vals = [d.output.vals(d.output.time <= cutoff_time,:);i3.vals(i3.time > cutoff_time,:) ];
		
		d.output = ra.TimeSeries(new_time,new_vals);
		d.prev = new_ts.burn(new_ts.tspan-10); % Retain only the last 5 seconds of the previous timeseries, should decrease memory load by a 

	end

end

function ts3 = append(ts1,ts2)
    % Append timeseries, assuming ts2 starts at the same time as ts1
    assert(max(abs(ts1.vals(end,:) - ts2.vals(1,:)))<1e-9,'ts1 end values do not match ts2 start values');
    tmp = ts2.remove_times(1); % Remove the duplicate time index
    tmp.time = tmp.time - tmp.time(1) + 1/tmp.fs + ts1.time(end); % Put it at the end of ts1
    ts3 = ra.TimeSeries([ts1.time;tmp.time],[ts1.vals;tmp.vals]);
end
		

