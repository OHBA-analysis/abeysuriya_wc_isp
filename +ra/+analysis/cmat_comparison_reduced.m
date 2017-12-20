function [model_cmat,band_frequency,state_corr] = cmat_comparison(wc,method,varargin)
	% Compare model connectivity to data - working at the moment only with alpha
	%
	% INPUTS
	% - wc - A simulated model run
	% - methods - One of 'aec','plv','pli'
	% - key-value pairs, see inputParser code
	% 	- aec_parent - handle to figure/container for plotting matrices
	% 	- corr_parent - handle to figure/container for plotting cross-frequency correlations
	% 	- orthongonalize - if true, raw time series will be orthogonalized prior to enveloping
	% 	- quick - if true, only the alpha-alpha correlation will be computed
	%
	% OUTPUTS
	% - model_cmat - cell array of connectivity matrices
	% - band_frequency - center frequencies for the bands (assuming +-2Hz width)
	% - state_corr - Cross frequency connectivity correlation matrix between model and data
	% - data - experimental connectivity matrices that were used (normally loaded from disk)

	arg = inputParser;
	arg.addParameter('orthogonalize',false);
	arg.addParameter('downsample',1); % Select downsampled frequency for AEC (does not apply to PLI or PLV)
	arg.addParameter('data',[]); % If empty, load data from file
	arg.parse(varargin{:});

	assert(ismember(method,{'aec','plv','pli'}));
	
	if isempty(arg.Results.data)
		data = load(fullfile(startup.get_rootdir,'romesh_nsys','data_files','data_fc'));
	else
		data = arg.Results.data;
	end
	
	switch method
		case 'aec'
			data_cmat = data.aecOrth; % AEC must be orthogonalized in data
		case 'plv'
			data_cmat = data.plvOrth; % PLV must be orthogonalized in data
		case 'pli'
			% PLI may or may not need orthogonalization - be consistent here
			if arg.Results.orthogonalize
				data_cmat = data.pliOrth;
			else
				data_cmat = data.pli;
			end
	end

	% Adjust frequencies if requested
	ts = wc.excitatory.copy;
	ts.cast('double'); % Make sure data type is correct

	if ts.tspan == 525
		fprintf(2,'Truncating due to longer than 500s due to 25s burn not being removed\n');
		ts.trim([15 10]);
	end

	connectivity_function = str2func(sprintf('ra.analysis.%s',method)); % Get the function to compute connectivity

	switch method 
		case 'aec'
			hts = ts.envelope('filter',data.band_frequency,'downsample',arg.Results.downsample,'orthogonalize',arg.Results.orthogonalize);
		case {'plv','pli'}
			[~,hts] = ts.envelope('filter',data.band_frequency,'orthogonalize',arg.Results.orthogonalize); % Note no downsampling for phase timeseries
	end

	model_cmat = connectivity_function(hts);
	band_frequency = data.band_frequency;

	% If we aren't using the DK parcellation, no point comparing to data
	if wc.excitatory.n_signals ~= size(data_cmat,1)
		model_cmat = nan(size(data_cmat,1));
		model_freqs = NaN;
		state_corr = NaN;
		return
	end

	% Compute correlation with data
	flt = find(triu(ones(size(model_cmat,1)),1)); % Upper triangle entries
	d = mean(data_cmat,3);
	state_corr = corr(d(flt),model_cmat(flt));
