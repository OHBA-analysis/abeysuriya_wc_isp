function mark_get_timeseries()
	% Starting with preprocessed but unfiltered sensor space MEEG
	% Filter it into a band (could be broadband)
	% And save the output file
	p = parcellation('/Users/romesh/oxford_postdoc/jonathan/data/nsys/bundle/stam_cortical.mat');
	subjects = {'3004','3006','3008','3013','3014','3015','3016','3018','3019','3020','3025','3029','3030','3031','3032','3033','3034','3035','3036','3037','3038','3039','3040','3044','3046','3047','3048','3049','3050','3051','3052','3053','3054','3055','3056','3057','3058','3060','3061','3062','3063','3064','3065','3068','3069','3070','3071','3072','3073','3074','3075','3076','3078','3080','3081'};
	parfor j = 1:length(subjects)
		D = spm_eeg_load(sprintf('/Users/romesh/oxford_postdoc/data/mark/beamformed/ffd%sopt_eo',subjects{j}));

		try
			D = D.montage('remove',3);
		end

		D.montage('switch',2);
		D = ROInets.get_node_tcs(D,p.parcelflag,'PCA');
		D.save()

	end

