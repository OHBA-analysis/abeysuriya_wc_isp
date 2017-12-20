function fig_9_cie

	d = load(fullfile('/Users/romesh/oxford_postdoc/romesh_nsys/sweeps','deco_rk4_adjustmean','analysis.mat'));
	d.inputs = d.inputs.simdata.inputs;
	r1 = load(fullfile('/Users/romesh/oxford_postdoc/romesh_nsys/sweeps','deco_rk4_adjustmean','results_part_1.mat'));
	[~,distance] = r1.out.wc.network.netmats;
	mean_distance = mean(distance(logical(triu(ones(size(distance)),1))));
	delay = 1000*mean_distance./d.inputs.velocity;
	coupling = d.inputs.coupling; 

	runs = ra.paper.individual_runs;
	runs.isp = load(runs.isp);
	hold on;
	scatter(runs.isp.out.wc.param.coupling(1),runs.isp.out.wc.mean_delay*1000,100,'ro');
	
	hfig = figure;
	subplot(1,3,1)
	imagesc(coupling, delay, d.outputs.cie_mean);
	cb = colorbar;
	xlabel('Coupling')
	ylabel('Mean delay (ms)')
	set(gca,'XTick',0.02:0.04:0.26)

	title('Homogeneous local inhibition')
	ylabel(cb,'$$c_{ie}$$','Interpreter','latex')
	axis square 

	subplot(1,3,2)
	imagesc(coupling, delay, d.outputs.alpha_synchrony);
	cb = colorbar;
	xlabel('Coupling')
	ylabel('Mean delay (ms)')
	set(gca,'XTick',0.02:0.04:0.26)

	title('Synchrony')
	ylabel(cb,'$$\overline{R}$$','Interpreter','latex')
	set(gca,'CLim',[0 1]);
	axis square 

	subplot(1,3,3)
	imagesc(coupling, delay, d.outputs.alpha_metastability);
	cb = colorbar;
	xlabel('Coupling')
	ylabel('Mean delay (ms)')
	set(gca,'XTick',0.02:0.04:0.26)

	title('Metastability')
	ylabel(cb,'$$\\std(R)$$','Interpreter','latex')
	set(gca,'CLim',[0 0.36]);
	axis square 

	set(gcf,'Position',[440   564   951   234]);

	% Get the run marker coordinates
	runs = ra.paper.individual_runs;
	runs.isp = load(runs.isp);
	subplot(1,3,1);hold on;scatter(runs.isp.out.wc.param.coupling(1),runs.isp.out.wc.mean_delay*1000,40,'ro');
	subplot(1,3,2);hold on;scatter(runs.isp.out.wc.param.coupling(1),runs.isp.out.wc.mean_delay*1000,40,'ro');
	subplot(1,3,3);hold on;scatter(runs.isp.out.wc.param.coupling(1),runs.isp.out.wc.mean_delay*1000,40,'ro');


	ra.plot.pfig(hfig,'fig_9b.pdf')





