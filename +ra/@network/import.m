function B = import(fname)
	% Load a bundle and convert it to a graph
	if ~exist(fname)
		fname = fullfile(getenv('NSYS_ROOT'),'../','Jonathan','data','nsys','bundle',fname);
	end
	
	d = load(fname);
	B = ra.network(d.conn,d.dist,d.short_names,d.long_names,d.cortical,d.roi_coords);
end
