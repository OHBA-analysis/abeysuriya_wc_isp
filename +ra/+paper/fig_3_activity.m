function fig_3_activity
	activity_timeseries
	ra.plot.pfig('fig_3_spectrum.pdf')
	close
	ra.plot.pfig('fig_3_timeseries.pdf')
	close

function activity_timeseries
	d = load(fullfile('/Users/romesh/oxford_postdoc/romesh_nsys/sweeps/deco_rk4_isp/analysis.mat'));

	r = ra.paper.individual_runs;
	r.isp = load(r.isp);

	par = parcellation('dk_cortical.nii.gz');
	idx = find(~cellfun(@isempty,strfind(par.labels,'Peri')));

	ts = r.isp.out.wc.excitatory;
	ts.cast('double');

	figure

	p(1) = plot(ts.select_signals(idx(1))) % Left pericalcirene

	[pks1,locs1] = findpeaks(ts.vals(:,idx(1)));
	[pks2,locs2] = findpeaks(ts.vals(:,idx(2)));

	Hen = ts.envelope('filter',[8 13],'orthogonalize',true);

	% Because the signal is orthogonalized, the power is potentially changed
	% Fit a scaling constant here
	scale_1 = polyfit(Hen.vals(100:end-100,idx(1)),interp1(ts.time(locs1),pks1,ts.time(100:end-100),'linear'),1);
	scale_2 = polyfit(Hen.vals(100:end-100,idx(2)),interp1(ts.time(locs2),pks2,ts.time(100:end-100),'linear'),1);

	hold on

	p(2) = plot(Hen.time,Hen.vals(:,idx(1))*scale_1(1) + scale_1(2),'r','LineWidth',2);
	p(3) = plot(Hen.time,Hen.vals(:,idx(2))*scale_2(1) + scale_2(2),'g','LineWidth',2);

	set(gca,'XLim',[275 275+30],'YLim',[0 1]);
	set(gca,'XTick',275:5:(275+30),'XTickLabel',[275:5:275+30]-275)
	set(gca,'YTick',0:0.2:1)
	set(gca,'Fontsize',16)
	ylabel('Normalised firing rate')
	uistack(p(1),'top');

	set(gcf,'Position',[   440   554   765   244]);

	% First, compute mean power spectrum
	P = [];
	for j = 1:ts.ns
		[P(:,j),f] = pwelch(detrend(ts.vals(:,j)),round(10*ts.fs),[],[],ts.fs);
	end
	P = mean(P,2);

	figure
	loglog(f,P)
	hold on
	fbar = @(x) plot([x x],[1e-11 1e-1],'r');
	arrayfun(fbar,[8 13]);
	delete(findobj(gcf,'Type','legend'))
	set(gca,'YLim',[1e-9,1e-1],'XLim',[1 100])
	set(gcf,'Position',[1081 553 290 242]);
	set(gca,'Fontsize',16)
	set(gca,'YTick',10.^[-8:2:-2])
	xlabel('Frequency (Hz)');
	ylabel('Power spectral density (s)')

