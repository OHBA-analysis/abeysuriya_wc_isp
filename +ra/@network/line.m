function B = line(n_units)
	% Return connection and distance matrices for a line of ROIs

	dmat = zeros(n_units);
	dmat(1,:) = 0:n_units-1;
	for j = 1:n_units
		dmat(j,:) = abs(dmat(1,:)-j+1);
	end

	cmat = 1./dmat;
	cmat(logical(eye(n_units)))=0;
	
	m = nsl.math.norm.amedian(dmat,true);
	%dmat = dmat*0.07; % Median full network distance
	dmat = dmat*0.0225/m; % Median subcortical distance

	B = ra.network(cmat,dmat);