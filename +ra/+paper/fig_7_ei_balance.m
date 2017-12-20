function ei_balance(run_name)
	if nargin < 1 || isempty(run_name) 
		run_name = 'deco_rk4_isp';
	end
	

	d = load(fullfile('/Users/romesh/oxford_postdoc/romesh_nsys/sweeps',run_name,'analysis.mat'));
	d.inputs = d.inputs.simdata.inputs;

	r = ra.paper.individual_runs;
	r1 = load(r.isp);
	[~,distance] = r1.out.wc.network.netmats;
	mean_distance = mean(distance(logical(triu(ones(size(distance)),1))));
	
	% X and Y axes
	delay = 1000*mean_distance./d.inputs.velocity;
	coupling = d.inputs.coupling; 

	figure

	subplot(1,2,1)
	imagesc(coupling, delay, sqrt(d.outputs.ei_rsquared));
	set(gca,'CLim',[0 1]);
	cb = colorbar
	xlabel('Coupling')
	ylabel('Mean delay (ms)')
	set(gca,'XTick',0.02:0.08:0.26)

	ylabel(cb,'Correlation (magnitude)')
	title('E/I balance')
	set(gca,'FontSize',14)
	axis square 

	subplot(1,2,2)
	scatter(r1.out.wc.network.strength.',r1.out.wc.cie,30,'bo','filled');
	xlabel('Network node strength')
	ylabel('c_{ie}')
	set(gca,'Fontsize',14)
	box on
	grid on
	set(gca,'XLim',[2 11],'YLim',[ -3.5   -2.5])

	str = r1.out.wc.network.strength.';
	p = polyfit(str,r1.out.wc.cie,1)
	cv = polyval(p,str);
	yresid = r1.out.wc.cie - cv;
	SSresid = sum(yresid.^2);
	SStotal = (length(cv)-1) * var(r1.out.wc.cie);
	ei_rsquared = 1 - SSresid/SStotal;

	hold on
	text(10.9,-2.51,sprintf('R^2=%.2f',ei_rsquared),'HorizontalAlignment','Right','VerticalAlignment','top','FontSize',14)
	
	plot([2 11],polyval(p,[2 11]),'r','LineWidth',1)

	%plot([2 11],ones(2,1).*mean(r1.out.wc.cie),'r--','LineWidth',1)

	set(gcf,'Position',[440   546   889   252]);


	subplot(1,2,1);hold on;scatter(0.14,9,40,'ro','filled');

	ra.plot.pfig('fig_7_ei_balance.pdf')