function varargout = pfig_montage(fname,h)
	% Save all open figures into a single png
	if nargin < 2 || isempty(h) 
		h = [];
	end
	
	if nargin < 1 || isempty(fname) 
		fname = 'default';
	end
	
	% Save all figures, and tile them
	if isempty(h)
		h = get(0,'Children');
	end
	
	fn = {};
	for j = 1:length(h)
		fn{j}= ra.plot.pfig(h(j),'%s_%d.png',fname,j);
	end

	output_dir = fileparts(fn{1})
	system(sprintf('/usr/local/bin/montage %s/%s_*.png  -geometry +1+1 %s/%s.png',output_dir,fname,output_dir,fname));

	for j = 1:length(fn)
		delete(fn{j});
	end

