function p = analyser(fname)
	% For convenience, fname is assumed to be a folder name within Romesh storage
	fname = fullfile(getenv('NSYS_ROOT'),'sweeps',fname,'results.mat');
	p = ra.sweep.analyse(fname);
