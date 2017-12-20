function connectivity_maps
	d = load('/Users/romesh/oxford_postdoc/mark_data/mark_final_connectivity_alpha.mat');

	% Which bands to show?
	% 2-6, 8-12, 20-24 as a limited set in primary text because these have very different data connectivity patterns
	% All of them in supplementary
	% These bands are - indexes 1,4,10

	ra.plot.tile_cmats({d.aecOrth},'titles',cellfun(@(x) sprintf('%d-%dHz',x),{d.band_frequency},'uniformoutput',false),'builtin_subplot',true)
	set(gcf,'Position',[440   613   868   185]);
	title('AEC')
	ra.plot.pfig('fig_5_data_aec.pdf')

	ra.plot.tile_cmats({d.plvOrth},'titles',cellfun(@(x) sprintf('%d-%dHz',x),{d.band_frequency},'uniformoutput',false),'builtin_subplot',true)
	set(gcf,'Position',[440   613   868   185]);
	title('PLV')
	ra.plot.pfig('fig_5_data_plv.pdf')

	ra.plot.tile_cmats({d.pli},'titles',cellfun(@(x) sprintf('%d-%dHz',x),{d.band_frequency},'uniformoutput',false),'builtin_subplot',true)
	set(gcf,'Position',[440   613   868   185]);
	title('PLI')
	ra.plot.pfig('fig_5_data_pli.pdf')

	runs = ra.paper.individual_runs;
	isp_run = load(runs.isp);

	model_aec = ra.analysis.cmat_comparison_reduced(isp_run.out.wc,'aec','orthogonalize',true);
	ra.plot.tile_cmats({model_aec},'titles',cellfun(@(x) sprintf('%d-%dHz',x),{d.band_frequency},'uniformoutput',false),'builtin_subplot',true)
	set(gcf,'Position',[440   613   868   185]);
	title('AEC')
	set(gca,'CLim',[-0.4    0.75])
	ra.plot.pfig('fig_5_model_aec.pdf')

	model_plv = ra.analysis.cmat_comparison_reduced(isp_run.out.wc,'plv','orthogonalize',true);
	ra.plot.tile_cmats({model_plv},'titles',cellfun(@(x) sprintf('%d-%dHz',x),{d.band_frequency},'uniformoutput',false),'builtin_subplot',true)
	set(gcf,'Position',[440   613   868   185]);
	title('PLV')
	set(gca,'CLim',[0 0.6]);
	ra.plot.pfig('fig_5_model_plv.pdf')


	model_pli = ra.analysis.cmat_comparison_reduced(isp_run.out.wc,'pli','orthogonalize',false);
	ra.plot.tile_cmats({model_pli},'titles',cellfun(@(x) sprintf('%d-%dHz',x),{d.band_frequency},'uniformoutput',false),'builtin_subplot',true)
	set(gcf,'Position',[440   613   868   185]);
	title('PLI')
	set(gca,'CLim',[0 1]);
	ra.plot.pfig('fig_5_model_pli.pdf')


	flt = find(triu(ones(size(model_aec,1)),1)); % Upper triangle entries
	data_aec = mean(d.aecOrth,3);
	data_plv = mean(d.plvOrth,3);
	data_pli = mean(d.pli,3);
	corr(data_aec(flt),model_aec(flt))
	corr(data_plv(flt),model_plv(flt))
	corr(data_pli(flt),model_pli(flt))