function fig_6_synchrony
	f1 = synchrony_summary('deco_rk4');
	ra.plot.pfig(f1(1),'fig_6_synchrony_noisp.pdf')
	ra.plot.pfig(f1(2),'fig_6_metastability_noisp.pdf')

	f2 = synchrony_summary('deco_rk4_isp');
	ra.plot.pfig(f2(1),'fig_6_synchrony_isp.pdf')
	ra.plot.pfig(f2(2),'fig_6_metastability_isp.pdf')

	d = load(fullfile(startup.get_rootdir,'romesh_nsys','data_files','data_synchrony.mat'))
	f3 = synchrony_summary('deco_rk4_isp');
	set(findobj(f3(1),'Type','image'),'CData',d.synchronyOrth)
	set(findobj(f3(2),'Type','image'),'CData',d.metastabilityOrth)
	ra.plot.pfig(f3(1),'fig_6_synchrony_isp_orth.pdf')
	ra.plot.pfig(f3(2),'fig_6_metastability_isp_orth.pdf')

	f2 = synchrony_summary('deco_rk4_isp_target_10_hd');
	ra.plot.pfig(f2(1),'fig_14a_synchrony_isp.pdf')
	ra.plot.pfig(f2(2),'fig_14b_metastability_isp.pdf')

	f2 = synchrony_summary('deco_rk4_isp_target_30_hd');
	ra.plot.pfig(f2(1),'fig_14c_synchrony_isp.pdf')
	ra.plot.pfig(f2(2),'fig_14d_metastability_isp.pdf')

function f = synchrony_summary(run_name)

	d = load(fullfile('/Users/romesh/oxford_postdoc/romesh_nsys/sweeps',run_name,'analysis.mat'));
	d.inputs = d.inputs.simdata.inputs;
	r1 = load(fullfile('/Users/romesh/oxford_postdoc/romesh_nsys/sweeps',run_name,'results_part_1.mat'));
	[~,distance] = r1.out.wc.network.netmats;
	mean_distance = mean(distance(logical(triu(ones(size(distance)),1))));
	
	% X and Y axes
	delay = 1000*mean_distance./d.inputs.velocity;
	coupling = d.inputs.coupling; 

	f(1) = figure;
	imagesc(coupling, delay, d.outputs.alpha_synchrony);
	cb = colorbar;
	xlabel('Coupling')
	ylabel('Mean Delay (ms)')
	set(gca,'XTick',0.02:0.04:0.26)
	title('Synchrony')
	ylabel(cb,'$$\overline{R}$$','Interpreter','latex')
	set(gca,'CLim',[0 1]);
	set(gca,'FontSize',20)
	axis square 

	f(2) = figure;
	imagesc(coupling, delay, d.outputs.alpha_metastability);
	max(d.outputs.alpha_metastability(:))
	min(d.outputs.alpha_metastability(:))
	cb = colorbar;
	xlabel('Coupling')
	ylabel('Mean delay (ms)')
	set(gca,'XTick',0.02:0.04:0.26)

	title('Metastability')
	ylabel(cb,'$$\\std(R)$$','Interpreter','latex')
	set(gca,'CLim',[0 0.36]);
	set(gca,'FontSize',20)
	axis square 

	if strcmp(run_name,'deco_rk4_isp')
		runs = ra.paper.individual_runs;
		runs.isp = load(runs.isp);
		figure(f(1));
		hold on;
		scatter(runs.isp.out.wc.param.coupling(1),runs.isp.out.wc.mean_delay*1000,100,'ro','filled');
		figure(f(2))
		hold on;
		scatter(runs.isp.out.wc.param.coupling(1),runs.isp.out.wc.mean_delay*1000,100,'ro','filled');
	end

	if strcmp(run_name,'deco_rk4_adjustmean') || strcmp(run_name,'deco_rk4_isp_target_10_hd') || strcmp(run_name,'deco_rk4_isp_target_30_hd')
		runs = ra.paper.individual_runs;
		runs.isp = load(runs.isp);
		figure(f(1));
		hold on;
		scatter(runs.isp.out.wc.param.coupling(1),runs.isp.out.wc.mean_delay*1000,100,'ro');
		figure(f(2))
		hold on;
		scatter(runs.isp.out.wc.param.coupling(1),runs.isp.out.wc.mean_delay*1000,100,'ro');
	end




