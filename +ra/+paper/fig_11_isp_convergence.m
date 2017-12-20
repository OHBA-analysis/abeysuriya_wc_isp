function fig_9_isp_convergence
	f = isp_convergence
	ra.plot.pfig(f(1),'fig_11_variability.pdf')
	ra.plot.pfig(f(2),'fig_11_cie_comparison.pdf')
	ra.plot.pfig(f(3),'fig_11_cie_timecourse.png')
	ra.plot.pfig(f(4),'fig_11_activity_timecourse.png')


function f = isp_convergence
	run_name = 'deco_rk4_isp';

	runs = ra.paper.individual_runs;
	high_synchrony = load(runs.high_synchrony);
	mid_synchrony = load(runs.isp);
	low_synchrony = load(runs.low_synchrony);

	d = load(fullfile('/Users/romesh/oxford_postdoc/romesh_nsys/sweeps',run_name,'analysis.mat'));
	d.inputs = d.inputs.simdata.inputs;
	r1 = load(fullfile('/Users/romesh/oxford_postdoc/romesh_nsys/sweeps',run_name,'results_part_1.mat'));
	[~,distance] = r1.out.wc.network.netmats;
	mean_distance = mean(distance(logical(triu(ones(size(distance)),1))));
	
	% X and Y axes
	delay = 1000*mean_distance./d.inputs.velocity;
	coupling = d.inputs.coupling; 

	f(1) = figure;
	imagesc(coupling, delay, d.outputs.mean_isp_variability);
	cb = colorbar;
	ylabel(cb,'c_{ie} variability')
	xlabel('Coupling')
	ylabel('Mean delay (ms)')
	set(gca,'XTick',0.02:0.04:0.26)

	axis square 
	set(gca,'Fontsize',20)

	t(1) = text(high_synchrony.out.wc.param.coupling(1),high_synchrony.out.wc.mean_delay*1000,'(1)','Color','r','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',20,'FontWeight','bold')
	t(2) = text(mid_synchrony.out.wc.param.coupling(1),mid_synchrony.out.wc.mean_delay*1000,'(2)','Color','r','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',20,'FontWeight','bold')
	t(3) = text(low_synchrony.out.wc.param.coupling(1),low_synchrony.out.wc.mean_delay*1000,'(3)','Color','r','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',20,'FontWeight','bold')

	f(3) = figure;
	ax(1) = subplot(1,3,1)
	high_synchrony.out.isp_ts.plot
	ylabel('c_{ie}')
	ra.plot.fig_letter('(f) - 1');

	ax(2) = subplot(1,3,2)
	mid_synchrony.out.isp_ts.plot
	ylabel('c_{ie}')
	ra.plot.fig_letter('(g) - 2');

	ax(3) = subplot(1,3,3)
	low_synchrony.out.isp_ts.plot
	ylabel('c_{ie}')
	ra.plot.fig_letter('(h) - 3');

	set(ax,'Fontsize',16)
	set(gcf,'Position',[157         313        1167         285])


	hs = high_synchrony.out.wc.iwe;
	ms = mid_synchrony.out.wc.iwe;
	ls = low_synchrony.out.wc.iwe;
	
	f(2) = figure('Position',[440   521   442   277])
	p(1) = plot([0 70],[0.15 0.15],'r');
	hold on
	p(2) = plot(1:68,sort(hs),'k--');
	p(3) = plot(1:68,sort(ms),'k');
	p(4) = plot(1:68,sort(ls),'k:');
	set(p,'LineWidth',2)
	legend('ISP target','1 - High synchrony','2 - Optimal parameters','3 - Low synchrony','location','NorthWest')
	set(gca,'FontSize',16)
	xlabel('ROI ')
	set(gca,'XTick',[1 15 34 50 68])
	set(gca,'XLim',[1 68])
	ylabel('Normalised firing rate')

	fs = high_synchrony.out.wc.excitatory.fs;
	f(4) = figure
	ax(1) = subplot(1,3,1)
	high_synchrony.out.wc.excitatory.select_times([1:(60*fs)]).plot
	ylabel('Normalised firing rate')
	set(gca,'YLim',[0 1])
	ra.plot.fig_letter('(c) - 1');

	ax(2) = subplot(1,3,2)
	mid_synchrony.out.wc.result.burn(120); % Of the 500 seconds, show a more interesting/typical segment
	mid_synchrony.out.wc.excitatory.select_times([1:(60*fs)]).plot
	ylabel('Normalised firing rate')
	set(gca,'YLim',[0 1])
	ra.plot.fig_letter('(d) - 2');

	ax(3) = subplot(1,3,3)
	low_synchrony.out.wc.excitatory.select_times([1:(60*fs)]).plot
	ylabel('Normalised firing rate')
	set(gca,'YLim',[0 1])
	ra.plot.fig_letter('(e) - 3');

	set(ax,'Fontsize',16)
	set(gcf,'Position',[157         313        1167         285])
