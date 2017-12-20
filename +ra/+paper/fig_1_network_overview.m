function network_overview

	% Plot the cortical surface
	p = parcellation('dk_cortical.nii.gz');
	net = ra.network.import('stam_cortical');

	[conn,dist] = net.netmats;

	p.plot_network(conn,0.9,60);
	set(gca,'CLim',[0 1]);
	colormap(bluewhitered(256))

	set(gca,'View', [0 90]);

	cb = findobj(gcf,'Type','ColorBar');
	set(gca,'FontSize',14)
	ylabel(cb,'Connection weight','FontSize',16);
	set(cb,'Position',[0.7670    0.2410    0.0357    0.5526])
	ra.plot.pfig('network_brain.pdf')

	figure
	imagesc(ra.analysis.reorder_matrix(conn))
	xlabel('ROI')
	ylabel('ROI')
	axis equal
	axis tight
	cb = colorbar
	ylabel(cb,'Connection weight','FontSize',22);
	set(gca,'FontSize',20)
	ra.plot.pfig('network_cmat.pdf')

	figure
	imagesc(ra.analysis.reorder_matrix(dist*100))
	xlabel('ROI')
	ylabel('ROI')
	axis equal
	axis tight
	cb = colorbar
	ylabel(cb,'Distance (cm)','FontSize',22);
	set(gca,'FontSize',20)
	ra.plot.pfig('network_dmat.pdf')


	% figure
	% imagesc(ra.analysis.reorder_matrix(d2/10))
	% xlabel('ROI')
	% ylabel('ROI')
	% axis equal
	% axis tight
	% cb = colorbar
	% ylabel(cb,'Distance (cm)');
	% set(gca,'FontSize',20)





	% c = p.roi_centers;
	% for j = 1:size(c,1)
	% 	d2(:,j) = sqrt(sum((c-c(j,:)).^2,2));
	% end
