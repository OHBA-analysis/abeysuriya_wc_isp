ra.paper.fig_4_fc_similarity
ax = findobj(0,'type','axes')
fo = findobj(0,'type','figure')
f = [];
for j = 1:length(ax)
f(j) = figure;
set(ax(j),'Parent',f(j));
end
close(fo)

set(ax,'FontSize',20,'Position',[0.2 0.2 0.6 0.6])
set(findobj(0,'type','scatter'),'SizeData',60)

for j = 1:length(ax)
	ra.plot.pfig(f(j),'p4_fc_%d.pdf',j)
	saveas(f(j),sprintf('~/Desktop/p4_fc_%d.fig',j))
end


runs = ra.paper.individual_runs;
high_synchrony = load(runs.high_synchrony);
mid_synchrony = load(runs.isp);
low_synchrony = load(runs.low_synchrony);
ax = [];

figure
high_synchrony.out.isp_ts.plot
title('High synchrony c_{ie}')
ylabel('c_{ie}')
set(findobj(gca,'type','line'),'LineWidth',1)
set(gca,'FontSize',22)
ra.plot.pfig('p6_a.png')

figure
mid_synchrony.out.isp_ts.plot
title('Optimal c_{ie}')
ylabel('c_{ie}')
set(findobj(gca,'type','line'),'LineWidth',1)
set(gca,'FontSize',22)
ra.plot.pfig('p6_b.png')


fs = high_synchrony.out.wc.excitatory.fs;
figure
high_synchrony.out.wc.excitatory.select_times([1:(60*fs)]).plot
title('High synchrony activity')
ylabel('Excitatory activity')
set(gca,'YLim',[0 1])
set(gca,'FontSize',22)
ra.plot.pfig('p6_c.png')

figure
mid_synchrony.out.wc.result.burn(120); % Of the 500 seconds, show a more interesting/typical segment
mid_synchrony.out.wc.excitatory.select_times([1:(60*fs)]).plot
title('Optimal activity')
ylabel('Excitatory activity')
set(gca,'YLim',[0 1])
set(gca,'FontSize',22)
ra.plot.pfig('p6_d.png')


