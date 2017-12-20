function fig_10_isp_target_convergence
	isp_nonconverge_new

	ra.plot.pfig('fig_12_timeseries_comparison.png')
	close
	ra.plot.pfig('fig_12_variability.pdf')

function isp_nonconverge_new

	d(1) = load(fullfile('/Users/romesh/oxford_postdoc/romesh_nsys/sweeps','deco_rk4_isp_target_10_hd','analysis.mat'));
	d(2) = load(fullfile('/Users/romesh/oxford_postdoc/romesh_nsys/sweeps','deco_rk4_isp','analysis.mat'));
	d(3) = load(fullfile('/Users/romesh/oxford_postdoc/romesh_nsys/sweeps','deco_rk4_isp_target_30_hd','analysis.mat'));
	% d(4) = load(fullfile('/Users/romesh/oxford_postdoc/romesh_nsys/sweeps','deco_rk4_isp_target_65','analysis.mat'));

	t = {'\rho = 0.1','\rho = 0.15','\rho = 0.3','\rho = 0.65'};

	runs = ra.paper.individual_runs;
	runs.high_isp_worst = load(runs.high_isp_worst);
	runs.isp = load(runs.isp);

	[~,distance] = runs.high_isp_worst.out.wc.network.netmats;
	mean_distance = mean(distance(logical(triu(ones(size(distance)),1))));

	figure('position',[137    85   817   558])
	set(gcf,'Position',[  137         417        1155         226])
	for j = 1:3
		subplot(1,3,j)
		delay = 1000*mean_distance./d(j).inputs.simdata.inputs.velocity;
		imagesc(d(j).inputs.simdata.inputs.coupling, delay, d(j).outputs.mean_isp_variability);
		cb = colorbar;
		ylabel(cb,'c_{ie} variability')
		xlabel('Coupling')
		ylabel('Mean delay (ms)')
		set(gca,'XTick',0.02:0.08:0.26)
		title(t{j})
		axis square 
		set(gca,'Fontsize',16)
		set(gca,'CLim',[0 0.11])

		if j == 2
			hold on;scatter(runs.isp.out.wc.param.coupling(1),runs.isp.out.wc.mean_delay*1000,40,'ro','filled');
		end

		if j == 3
			hold on;scatter(runs.high_isp_worst.out.wc.param.coupling(1),runs.high_isp_worst.out.wc.mean_delay*1000,40,'ro','filled');
		end
	end


	figure('position',[257   577   915   219])
	dr = load('/Users/romesh/oxford_postdoc/romesh_nsys/sweeps/isp_nonconverge_detailed_example/result_run.mat')
	subplot(1,3,1)
	plot(runs.high_isp_worst.out.isp_ts)
	ylabel('c_{ie}')
	set(gca,'FontSize',16)
	subplot(1,3,2)
	plot(dr.isp_ts)
	ylabel('c_{ie}')
	set(gca,'FontSize',16)
	set(gca,'XLim',[50 300],'XTick',50:50:300,'XTickLabel',(50:50:300)-50)
	subplot(1,3,3)
	plot(dr.wc.excitatory)
	ylabel('Normalised firing rate')
	set(gca,'FontSize',16)
	set(gca,'XLim',[50 300],'XTick',50:50:300,'XTickLabel',(50:50:300)-50)

	return

	figure
	x = dr.wc.excitatory;
	xi = dr.wc.inhibitory;

	y = dr.isp_ts;
	window = 0.25;

	clear xs
	for j = 1:68
		xs(:,j) = smooth(x.vals(:,j),window*x.fs);
	end
	scatter(mean(y.vals,2),mean(xs,2))
	xlabel('Mean c_{ie}')
	ylabel('Mean excitatory activity')
	box on
	set(gca,'FontSize',16)


	% Compute smoothed IWE
	winlen = 600;
	idx = 1:(size(x.vals,1)-winlen);
	output = zeros(length(idx),size(x.vals,2));
	tic
	for j = 1:length(idx)
		a = idx(j):idx(j)+winlen;
		output(j,1) = sum(x.vals(a,1).*y.vals(a,1))./sum(y.vals(a,1));
	end
	toc
	% scatter(mean(y.vals(idx+winlen/2,:),2),mean(output,2))

	scatter(mean(y.vals(idx+winlen/2,1),2),mean(output(:,1),2))

	return
	% 	figure
	% y2 = mean(y.vals(truncate*y.fs:end-truncate*y.fs,:),2);
	% for j = 1:68
	% 	xs2(:,j) = smooth(x.vals(truncate*x.fs:end-truncate*x.fs,j),150);
	% 	[a{j},b{j}] = findpeaks(xs2(:,j));
	% 	hold on
	% 	scatter(y2(b{j}),a{j})
	% end


	% for j = 1:





