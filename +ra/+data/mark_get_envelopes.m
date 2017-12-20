function mark_get_envelopes(dfilt)
	% Read the beamformed MEEG objects
	% Apply filtering, parcellation, enveloping, and downsampling

	if nargin < 1 || isempty(dfilt) 
		dfilt = false; % Filter D object prior to doing PCA
	end
	
	band_names = {'alpha'};
	band_frequencies = {[8 13]};

	for j = 1:length(band_names)
		fprintf('Processing band %s\n',band_names{j});
		mark_process_band(dfilt,band_names{j},band_frequencies{j});
	end
end

function mark_process_band(dfilt,band_name,band_frequency)
	% Define control subjects
	p = parcellation('dk_cortical.nii.gz');
	subjects = {'3004','3006','3008','3013','3014','3015','3016','3018','3019','3020','3025','3029','3030','3031','3032','3033','3034','3035','3036','3037','3038','3039','3040','3044','3046','3047','3048','3049','3050','3051','3052','3053','3054','3055','3056','3057','3058','3060','3061','3062','3063','3064','3065','3068','3069','3070','3071','3072','3073','3074','3075','3076','3078','3080','3081'};
	%subjects = subjects(1);

	aec = nan(p.n_parcels,p.n_parcels,length(subjects));
	aecOrth = aec;
	plv = aec;
	plvOrth = aec;
	pli = aec;
	pliOrth = aec;

	for j = 1:length(subjects)
		
		fprintf('Subject %d of %d\n',j,length(subjects));

		% Load the beamformed timeseries
		D = spm_eeg_load(sprintf('/Users/romesh/oxford_postdoc/data/mark/beamformed/ffd%sopt_eo',subjects{j}));

		if dfilt
			fprintf('Using MEEG filtering\n')
		
			% Filter it into this band
			working_file = tempname(pwd);
			D = D.montage('switch',0);
			D = D.copy(working_file);
			osl_spmfun(@spm_eeg_filter,struct('D',working_file,'band','bandpass','dir','twopass','freq',band_frequency));
			D = spm_eeg_load(working_file);

			% Switch to the beamformed montage
			D = D.montage('switch',2);

			% Get parcel timecourses in this frequency band
			D = ROInets.get_node_tcs(D,p.parcelflag,'PCA');

			% Get node timeseries
			ts = ra.TimeSeries(D.time,D(:,:,:).');
			clean = all(D.badsamples(:,:,:)==0,1);

			D.delete();

		else
			fprintf('Using ts filtering\n')

			D = spm_eeg_load(sprintf('/Users/romesh/oxford_postdoc/data/mark/beamformed/ffd%sopt_eo',subjects{j}));
			D = D.montage('switch',3);
			ts = ra.TimeSeries(D.time,D(:,:,:).');
			clean = all(D.badsamples(:,:,:)==0,1);
			ts.filter(band_frequency); % Filter it
		end

		% Make the envelopes
		[Hen(j),HPhase(j)] = ts.envelope(); % Don't rescale, don't downsample
		plv(:,:,j) = ra.analysis.plv(HPhase(j).select_times(clean));
		pli(:,:,j) = ra.analysis.pli(HPhase(j).select_times(clean));

		ts.orthogonalize; % Apply orthogonalization

		% Get the PLV metrics again
		[HenOrth(j),HPhaseOrth(j)] = ts.envelope(); % Don't rescale, don't downsample
		plvOrth(:,:,j) = ra.analysis.plv(HPhaseOrth(j).select_times(clean));
		pliOrth(:,:,j) = ra.analysis.pli(HPhaseOrth(j).select_times(clean));

		% Now assign the downsampled outputs
		Hen1(j) = ds_clean(Hen(j),1,clean);
		HenOrth1(j) = ds_clean(HenOrth(j),1,clean);

		aec(:,:,j) = ra.analysis.aec(Hen1(j));
		aecOrth(:,:,j) = ra.analysis.aec(HenOrth1(j));

    end
    
	save(sprintf('mark_data/mark_final_%s',band_name),'dfilt','band_frequency','aec','plv','pli','aecOrth','plvOrth','pliOrth','Hen','Hen1','HenOrth','HenOrth1','HPhase','HPhaseOrth','subjects','-v7.3');
	save(sprintf('mark_data/mark_final_connectivity_%s',band_name),'dfilt','band_frequency','aec','plv','pli','aecOrth','plvOrth','pliOrth','subjects','-v7.3');

end

function Hen2 = ds_clean(Hen,fs,clean)
	% Now assign the downsampled outputs
	Hen2 = Hen.resample(fs);
	mapping = ra.map_times(Hen.time,Hen2.time);
	Hen2.vals(mapping(~clean),:) = NaN;
end
