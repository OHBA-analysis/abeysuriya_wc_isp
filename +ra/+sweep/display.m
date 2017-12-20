function display(run_name)

	d = load(fullfile('/Users/romesh/oxford_postdoc/romesh_nsys/sweeps',run_name,'results.mat'));
	d2 = load(fullfile('/Users/romesh/oxford_postdoc/romesh_nsys/sweeps',run_name,'analysis.mat'));
	d.outputs = d2.outputs;
	r1 = load(fullfile('/Users/romesh/oxford_postdoc/romesh_nsys/sweeps',run_name,'results_part_1.mat'));


	if size(get(0, 'MonitorPositions'),1) > 1 
		figure('Position',[ 1666 -111 1481 1258 ],'Name',run_name)
	else
		figure('Position',[ 48 60 1280 680 ],'Name',run_name)
	end

	[~,distance] = r1.out.wc.network.netmats;
	mean_distance = mean(distance(logical(triu(ones(size(distance)),1))));

	if size(d.task_ids,3) > 1 % If is is an ISP sweep
		xvals = squeeze(d.inputs.isp_target(1,:,:));
		yvals = 1000*mean_distance./d.inputs.v(:,:,1);
		xlabelstring = 'ISP target';
		ylabelstring = 'Delay (ms)';
	else
		d.inputs.delay = 1000*mean_distance./d.inputs.velocity;
		xvals = d.inputs.coupling;
		yvals = d.inputs.delay;
		xlabelstring = 'Coupling';
		ylabelstring = 'Delay (ms)';
	end


	function standard_plot(datastr,titlestr,clim)
		% Take in the field name for struct d.outputs, or a matrix of actual values
		if nargin < 3 || isempty(clim) 
			clim = [];
		end
		
		ax(counter) = ra.plot.subplot(4,4,counter);
		marker(counter) = NaN;

		try
			if ischar(datastr)
				dx = d.outputs.(datastr);
			else
				dx = datastr;
			end

			imagesc(xvals,yvals,squeeze(dx),'HitTest','off','parent',ax(counter));
			hold(ax(counter),'on')
			marker(counter) = scatter(ax(counter),NaN,NaN,50,'ro','filled','hittest','off');
			xlabel(ax(counter),xlabelstring)
			ylabel(ax(counter),ylabelstring)
			title(ax(counter),titlestr)
			set(ax(counter),'XLim',[min(xvals) max(xvals)],'YLim',[min(yvals) max(yvals)])
			axis(ax(counter),'square');
			axis(ax(counter),'tight');
			colorbar('peer',ax(counter));
			set(ax(counter),'FontSize',12)
			%colormap jet
			if ~isempty(clim)
				set(ax(counter),'CLim',clim);
			end
		catch ME
			axis off
			fprintf(2,'%s\n',ME.message);
		end

		counter = counter + 1;

	end

	ax = [];
	marker = [];
	counter = 1;

	% First row - metastability metrics
	standard_plot('f_adjust_factor','Frequency scale factor');
	%standard_plot('power_bw','Power bandwidth (10% max)');
	standard_plot('power_std','Power variability (std)');
	standard_plot('f_max','Peak frequency (avg)');

	% Second row, asymmetry and oscillation
	%standard_plot('alpha_conn_aec','Alpha AEC '); % If peak frequency is matched, how good is the fit to alpha AEC?
	standard_plot((d.outputs.('alpha_conn_aec')-0.5794)./0.1829,'Alpha AEC zscore'); % If peak frequency is matched, how good is the fit to alpha AEC?
	standard_plot('alpha_conn_aec','Alpha AEC '); % If peak frequency is matched, how good is the fit to alpha AEC?

	standard_plot('alpha_conn_no_orth_aec','No orth Alpha AEC '); % If peak frequency is matched, how good is the fit to alpha AEC?
	standard_plot('alpha_conn_plv','Alpha PLV'); % If peak frequency is matched, how good is the fit to alpha AEC?
	standard_plot('alpha_conn_no_orth_plv','No orth Alpha PLV'); % If peak frequency is matched, how good is the fit to alpha AEC?

	standard_plot('alpha_conn_pli','Alpha PLI'); % If peak frequency is matched, how good is the fit to alpha AEC?
	standard_plot('alpha_conn_no_orth_pli','No orth Alpha PLI'); % If peak frequency is matched, how good is the fit to alpha AEC?
	standard_plot('ei_rsquared','EI R-squared'); % If peak frequency is matched, how good is the fit to alpha AEC?
	standard_plot('iwe_plastic_discrepancy','IWE target discrepancy');


	% % Third row, plasticity
	% standard_plot('max_isp_variability','ISP max variability');
	% standard_plot('iwe_discrepancy','IWE discrepancy');
	% standard_plot('excitatory_amp_discrepancy','Excitatory amplitude discrepancy');

	% % Fourth row, synchrony measures
	standard_plot('alpha_synchrony','Alpha synchrony',[0 1]);
	standard_plot('alpha_metastability','Alpha metastability');
	standard_plot('max_isp_variability','Max ISP variability');
	standard_plot('mean_isp_variability','Mean ISP variability');
	
	set(ax,'ButtonDownFcn',@(a,b,c) wrap(a,b,xvals,yvals,d,marker,run_name))

end

function wrap(a,b,cv,vv,d,marker,run_name)
	% Coordinates user clicked on
	clicked = b.IntersectionPoint(1:2);

	[~,cidx] = min(abs(cv-clicked(1)));
	[~,vidx] = min(abs(vv-clicked(2)));

	for j = 1:length(marker)
		try
			set(marker(j),'XData',cv(cidx),'YData',vv(vidx));
		end
	end
	drawnow
	
	idx = sub2ind([length(vv) length(cv)],vidx,cidx);
	fprintf('Run file: results_part_%d.mat\n',idx)
	dat = load(fullfile('/Users/romesh/oxford_postdoc/romesh_nsys/sweeps',run_name,sprintf('results_part_%d.mat',idx)));

	dat.out.wc.result.cast('double'); % Make results double automatically

	assignin('base','wc_selected',dat.out.wc)

	h = dat.out.wc.analyse;

	if isfield(dat.out,'isp_ts')
		userdata = get(h,'UserData');
		a = uitab(userdata.tg,'Title','ISP');
		b = uicontainer(a);
		ax = axes(b);
		plot(dat.out.isp_ts,'Parent',ax)
		assignin('base','wc_selected_isp_ts',dat.out.isp_ts);
	end

	if isfield(d.outputs,'alpha_conn_aec_matrix')
		userdata = get(h,'UserData');
		a = uitab(userdata.tg,'Title','AEC');
		b = uicontainer(a);
		ax = subplot(1,2,1,'Parent',b);
		imagesc(ra.analysis.reorder_matrix(d.outputs.alpha_conn_aec_matrix{idx}));
		title(ax,'Orthogonalized')
		colorbar('peer',ax)

		ax = subplot(1,2,2,'Parent',b);
		imagesc(ra.analysis.reorder_matrix(d.outputs.alpha_conn_no_orth_aec_matrix{idx}));
		title(ax,'No orthogonalization')
		colorbar('peer',ax)


		a = uitab(userdata.tg,'Title','PLV');
		b = uicontainer(a);
		ax = subplot(1,2,1,'Parent',b);
		imagesc(ra.analysis.reorder_matrix(d.outputs.alpha_conn_plv_matrix{idx}));
		title(ax,'Orthogonalized')
		colorbar('peer',ax)

		ax = subplot(1,2,2,'Parent',b);
		imagesc(ra.analysis.reorder_matrix(d.outputs.alpha_conn_no_orth_plv_matrix{idx}));
		title(ax,'No orthogonalization')
		colorbar('peer',ax)


		a = uitab(userdata.tg,'Title','PLI');
		b = uicontainer(a);

		ax = subplot(1,2,1,'Parent',b);
		imagesc(ra.analysis.reorder_matrix(d.outputs.alpha_conn_pli_matrix{idx}));
		title(ax,'Orthogonalized')
		colorbar('peer',ax)

		ax = subplot(1,2,2,'Parent',b);
		imagesc(ra.analysis.reorder_matrix(d.outputs.alpha_conn_no_orth_pli_matrix{idx}));
		title(ax,'No orthogonalization')
		colorbar('peer',ax)


		a = uitab(userdata.tg,'Title','Raw');
		b = uicontainer(a);

		ax = subplot(1,2,1,'Parent',b);
		imagesc(ra.analysis.reorder_matrix(d.outputs.raw_correlation{idx}));
		title(ax,'Raw correlation')
		colorbar('peer',ax)

		ax = subplot(1,2,2,'Parent',b);
		imagesc(diag(nan(68,1))+ra.analysis.reorder_matrix(d.outputs.raw_correlation_alpha{idx}));
		title(ax,'Raw correlation, alpha band')
		colorbar('peer',ax)

	end
end





