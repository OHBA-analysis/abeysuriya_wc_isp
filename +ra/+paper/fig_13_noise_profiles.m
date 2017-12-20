function fig_13_noise_profiles

	noises = {'noise_0x','noise_5x','noise_10x','noise_20x'};

	runs = ra.paper.individual_runs;


	for j = 1:length(noises)
		r = load(runs.(noises{j}));
		r.out.wc.units(1).E.P
		r.out.wc.units(1).I.P

		[model_aec,band_frequency] = ra.analysis.cmat_comparison_reduced(r.out.wc,'aec','orthogonalize',true);
		ra.plot.tile_cmats({model_aec},'titles',cellfun(@(x) sprintf('%d-%dHz',x),{band_frequency},'uniformoutput',false),'builtin_subplot',true)
		set(gca,'CLim',[-0.4    0.75])
		set(gcf,'Position',[440   613   868   185]);
		title('')
		xlabel('');
		ylabel('');

		if j > 1
			set(gca,'YTick',[]);
		end
		set(gca,'XTick',[]);

		if j < length(noises)
			delete(findobj(gcf,'type','colorbar'))
		end

		ra.plot.pfig('%s_aec.pdf',noises{j})

		model_plv = ra.analysis.cmat_comparison_reduced(r.out.wc,'plv','orthogonalize',true);
		ra.plot.tile_cmats({model_plv},'titles',cellfun(@(x) sprintf('%d-%dHz',x),{band_frequency},'uniformoutput',false),'builtin_subplot',true)
		set(gca,'CLim',[0 0.6]);
		set(gcf,'Position',[440   613   868   185]);
		title('')
		xlabel('');
		ylabel('');

		if j > 1
			set(gca,'YTick',[]);
		end
		set(gca,'XTick',[]);

		if j < length(noises)
			delete(findobj(gcf,'type','colorbar'))
		end

		ra.plot.pfig('%s_plv.pdf',noises{j})

		model_pli = ra.analysis.cmat_comparison_reduced(r.out.wc,'pli','orthogonalize',false);
		ra.plot.tile_cmats({model_pli},'titles',cellfun(@(x) sprintf('%d-%dHz',x),{band_frequency},'uniformoutput',false),'builtin_subplot',true)
		set(gcf,'Position',[440   613   868   185]);
		set(gca,'CLim',[0 1]);
		title('')
		xlabel('');
		ylabel('');		

		if j > 1
			set(gca,'YTick',[]);
		end

		if j < length(noises)
			delete(findobj(gcf,'type','colorbar'))
		end
		
		ra.plot.pfig('%s_pli.pdf',noises{j})

	end