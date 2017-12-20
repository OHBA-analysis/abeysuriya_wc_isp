function set_noise(wc,mu,sig)
	if nargin < 3 || isempty(sig) 
		sig = [];
	end
	
	if nargin < 2 || isempty(mu) 
		mu = [];
	end
	
	if ~isempty(mu) && length(mu) == 1
		mu = [mu mu];
	end

	if ~isempty(sig) && length(sig) == 1
		sig = [sig sig];
	end

	for j = 1:wc.n_units
		if ~isempty(mu)
			wc.units(j).E.P.mu = mu(1);
			wc.units(j).I.P.mu = mu(2);
		end
		if ~isempty(sig)
			wc.units(j).E.P.sigma = sig(1);
			wc.units(j).I.P.sigma = sig(2);
		end

		wc.units(j).E.P.name = 'gaussian';
		assert(isfield(wc.units(j).E.P,'mu'))
		assert(isfield(wc.units(j).E.P,'sigma'))
		wc.units(j).I.P.name = 'gaussian';
		assert(isfield(wc.units(j).I.P,'mu'))
		assert(isfield(wc.units(j).I.P,'sigma'))
	end

