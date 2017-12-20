function p = deco_rk4()
	% Set up the unit
	u = ra.unit.deco_10;
	noise_mu = [0.31 0];
	noise_sigma = [1e-2 1e-2];
	u.E.P = struct('mu',noise_mu(1),'sigma',noise_sigma(1),'name','gaussian');
	u.I.P = struct('mu',noise_mu(2),'sigma',noise_sigma(2),'name','gaussian');

	% Set up the network, global coupling, and delays
	net = ra.network.import('stam_cortical');
	coupling = 0.02:0.01:0.26;
	[c,d] = net.netmats;
	mean_distance = mean(d(logical(triu(ones(size(d)),1))));
	min_distance = min(d(logical(triu(ones(size(d)),1)))); % The maximum allowed velocity is min_distance/tstep. The smallest delay is then  mean_distance/(min_distance/tstep) 
	% For 5e-5, cortical has 
	mean_delays = 1e-3*(0:1:20);
	velocity = mean_distance./mean_delays;

	p = ra.sweep.simulate(ra.network.import('stam_cortical'),u,velocity,coupling);
