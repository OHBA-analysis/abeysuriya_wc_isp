function [synchrony,metastability] = synchrony(broadband_ts,freqs,orthogonalize,clean)
	% Broadband TS goes in because it needs to be filtered into each band being analyzed
	% varargin gets passed to envelope

	if nargin < 4 || isempty(clean) 
		clean = logical(ones(broadband_ts.n_times,1));
	end
	
	if nargin < 3 || isempty(orthogonalize) 
		orthogonalize = true;
	end

	if nargin < 2 || isempty(freqs) 
		[~,freqs] = ra.data.adam_bands(2);
	end
	
	synchrony = nan(size(freqs));
	metastability = nan(size(freqs));
	for j = 1:length(freqs)
		[~,Ph] = broadband_ts.envelope('filter',freqs{j},'orthogonalize',orthogonalize); % Don't downsample
		order_ts = abs(mean(exp(1i*Ph.vals(clean,:)),2));
		synchrony(j) = mean(order_ts);
		metastability(j) = std(order_ts);
	end
