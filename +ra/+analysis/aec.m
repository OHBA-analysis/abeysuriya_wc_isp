function c = aec(Hen)
	% Take in a envelope timeseries, compute correlation (zero diagonal)
	% If output isn't being captured, make a plot
	% If multiple timeseries are provided, average AEC over them
	
	aec = nan(Hen(1).n_signals,Hen(1).n_signals,length(Hen));

	if Hen(1).fs > 40
		fprintf('Frequency higher than 40Hz, did you envelope and downsample...?\n');
	end

	for j = 1:length(Hen)
		clean = all(isfinite(Hen(j).vals),2);
		aec(:,:,j) = corr(Hen(j).vals(clean,:)).*~eye(Hen(j).n_signals);
	end

	c = mean(aec,3);

	if nargout == 0
	    figure
	    imagesc(c.*~eye(size(c)));
	    ra.plot.format_corr;
	end
end
