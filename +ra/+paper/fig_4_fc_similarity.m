function fig_4_fc_similarity

	sweeps = {'deco_rk4','deco_rk4_isp','deco_rk4_adjustmean','deco_rk4_isp_target_10_hd','deco_rk4_isp_target_30_hd'};

	for j = 1:length(sweeps)
		h(j) = sweep_summary(sweeps{j});
	end

	for j = 1:3 % For each modality

		for k = 1:length(sweeps) % For each sweep
			figure(h(k));
			subplot(1,3,j);
			a(k) = min(get(gca,'CLim'));
			b(k) = max(get(gca,'CLim'));
		end

		for k = 1:length(sweeps)
			figure(h(k));
			subplot(1,3,j);
			set(gca,'CLim',[min(a) max(b)]);
		end

	end
		
	ra.plot.pfig(h(1),'fig_4a.pdf')
	ra.plot.pfig(h(2),'fig_4b.pdf')
	ra.plot.pfig(h(3),'fig_9a.pdf')
	ra.plot.pfig(h(4),'fig_8a.pdf')
	ra.plot.pfig(h(5),'fig_8b.pdf')

	sweep_summary('deco_rk4_isp',true)
	ra.plot.pfig('fig_4c.pdf')


function hfig = sweep_summary(run_name,do_zscore)
	if nargin < 2 || isempty(do_zscore) 
		do_zscore = false;
	else
		a = ra.data.mark_zscores;
	end
	
	d = load(fullfile('/Users/romesh/oxford_postdoc/romesh_nsys/sweeps',run_name,'analysis.mat'));
	d.inputs = d.inputs.simdata.inputs;
	r1 = load(fullfile('/Users/romesh/oxford_postdoc/romesh_nsys/sweeps',run_name,'results_part_1.mat'));
	[~,distance] = r1.out.wc.network.netmats;
	mean_distance = mean(distance(logical(triu(ones(size(distance)),1))));
	
	% X and Y axes
	delay = 1000*mean_distance./d.inputs.velocity;
	coupling = d.inputs.coupling; 

	hfig = figure;

	subplot(1,3,1)
	if do_zscore
		cx = (d.outputs.alpha_conn_aec-a.aecOrth(1))./a.aecOrth(2);
		imagesc(coupling, delay,cx);
		title('AEC')
		cb=colorbar;
		ylabel(cb,'Z')
		% set(cb,'Visible','off')
    else
		imagesc(coupling, delay, d.outputs.alpha_conn_aec);
		title('AEC')
		cb=colorbar;
		ylabel(cb,'Correlation')
	end

	xlabel('Coupling')
	ylabel('Mean delay (ms)')
	set(gca,'XTick',0.02:0.04:0.26)
	axis square 

	subplot(1,3,2)
	if do_zscore
		cx = (d.outputs.alpha_conn_plv-a.plvOrth(1))./a.plvOrth(2);
		imagesc(coupling, delay, cx);
		title('PLV')
		cb=colorbar;
		ylabel(cb,'Z')
		% set(cb,'Visible','off')
	else
		imagesc(coupling, delay, d.outputs.alpha_conn_plv);
		title('PLV')
		cb=colorbar;
		ylabel(cb,'Correlation')
	end
	xlabel('Coupling')
	ylabel('Mean delay (ms)')
	set(gca,'XTick',0.02:0.04:0.26)

	axis square 

	subplot(1,3,3)
	if do_zscore
		cx = (d.outputs.alpha_conn_no_orth_pli-a.pli(1))./a.pli(2);
		imagesc(coupling, delay, cx);
		title('PLI')
		cb=colorbar;
		ylabel(cb,'Z')
		% set(cb,'Visible','off')
	else
		imagesc(coupling, delay, d.outputs.alpha_conn_no_orth_pli);
		title('PLI')
		cb=colorbar;
		ylabel(cb,'Correlation')
	end
	xlabel('Coupling')
	ylabel('Mean delay (ms)')
	set(gca,'XTick',0.02:0.04:0.26)

	axis square 

	set(gcf,'Position',[440   564   951   234]);

	%return
	
	% Get the run marker coordinates
	runs = ra.paper.individual_runs;

	if strcmp(run_name,'deco_rk4_isp')
		runs.isp = load(runs.isp);
		subplot(1,3,1);hold on;scatter(runs.isp.out.wc.param.coupling(1),runs.isp.out.wc.mean_delay*1000,40,'ro','filled');
		subplot(1,3,2);hold on;scatter(runs.isp.out.wc.param.coupling(1),runs.isp.out.wc.mean_delay*1000,40,'ro','filled');
		subplot(1,3,3);hold on;scatter(runs.isp.out.wc.param.coupling(1),runs.isp.out.wc.mean_delay*1000,40,'ro','filled');
	elseif ~strcmp(run_name,'deco_rk4')
		runs.isp = load(runs.isp);
		subplot(1,3,1);hold on;scatter(runs.isp.out.wc.param.coupling(1),runs.isp.out.wc.mean_delay*1000,40,'ro');
		subplot(1,3,2);hold on;scatter(runs.isp.out.wc.param.coupling(1),runs.isp.out.wc.mean_delay*1000,40,'ro');
		subplot(1,3,3);hold on;scatter(runs.isp.out.wc.param.coupling(1),runs.isp.out.wc.mean_delay*1000,40,'ro');
	end

