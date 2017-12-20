function prepare_si

	% prepare parcellation file

	% prepare connectivity and delay matrices

	prefix = '/Users/romesh/Dropbox/work/oxford/paper-isp/supplementary';

	% Save the parcellation
	net = ra.network.import('stam_cortical');

	[cx,d]=net.netmats;
	namesx = net.nodes{:,'Longname'};

	[c,total_order] = ra.analysis.reorder_matrix(cx);
	[d] = ra.analysis.reorder_matrix(d);
	names = namesx(total_order);

	dlmwrite(fullfile(prefix,'connectivity.txt'),c,'precision','%.9f');
	dlmwrite(fullfile(prefix,'distance.txt'),d,'precision','%.9f');

	writetable(table(names),fullfile(prefix,'parcel_names.txt'),'WriteVariableNames',false);

	p = parcellation('dk_cortical.nii.gz');
	m = p.to_vol(p.value_vector);
	p.savenii(m,fullfile(prefix,'dk_parcellation.nii.gz'));

	p.weight_mask = p.weight_mask(:,:,:,total_order);
	m = p.to_vol(p.value_vector);
	p.savenii(m,fullfile(prefix,'dk_parcellation2.nii.gz'));

	