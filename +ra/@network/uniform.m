function B = uniform(n_units)
	% Return network with complete uniform connectivity
	dmat = ones(n_units);
	cmat = ones(n_units);
	cmat(logical(eye(n_units)))=0;
	B = ra.network(cmat,dmat);