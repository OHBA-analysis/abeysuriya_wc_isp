function frequency_comparison
	runs = ra.paper.individual_runs;
	isp_run = load(runs.isp);
	wc = isp_run.out.wc;
	wc.result.cast('double');
	ts = wc.excitatory;
	ts.trim([15 10]);
	%ts = ra.scale_frequency(ts,10,[0.5 2]); % Rescale frequency

	bands = {[8 13],[18 23],[28 33]}; % Traditional

	% Fully data-driven rather than starting with alpha
	% bc = 9.5704.*[1 2 3 4];
	% for j = 1:length(bc)
	% 	bands{j} = bc(j)+[-2 2];
	% end

	% ALSO spectrum
	% First, compute mean power spectrum
	P = [];
	for j = 1:ts.ns
		[P(:,j),f] = pwelch(detrend(ts.vals(:,j)),round(10*ts.fs),[],[],ts.fs);
	end
	P = mean(P,2);

	sf = figure
	loglog(f,P)
	af = gca;
	hold on
	fbar = @(x) plot([x x],[1e-11 1e-1],'r','parent',af);
	delete(findobj(gcf,'Type','legend'))
	set(gca,'YLim',[1e-8,0.5e-1],'XLim',[1 100])
	set(gcf,'Position',[   326   420   438   263]);
	set(gca,'Fontsize',16)
	set(gca,'YTick',10.^[-8:2:-2])
	xlabel('Frequency (Hz)');
	ylabel('Power spectral density (s)')

	for j = 1:length(bands)
		arrayfun(fbar,bands{j});
		Hen = ts.envelope('filter',bands{j},'downsample',1,'orthogonalize',true);
		ra.plot.tile_cmats({ra.analysis.aec(Hen)},'titles',{sprintf('%d-%dHz',bands{j}(1),bands{j}(2))},'builtin_subplot',true)
		set(gcf,'Position',[440   613   868   185]);
		ax(j) = gca;
		fh(j) = gcf;
	end

	cl(1) = min(arrayfun(@(x) min(get(x,'CLim')),ax));
	cl(2) = max(arrayfun(@(x) max(get(x,'CLim')),ax));
	%set(ax,'CLim',cl);

	for j = 1:length(bands)
		ra.plot.pfig(fh(j),'fig_10_%d.pdf',j)
	end


	ra.plot.pfig(sf,'fig_10_spec.pdf');

