function random_walk_driver
	% This function computes the expected synchrony based on a 2D random walk
	% with uniformly distributed direction and unit distance at each step.
	% The theoretical results come from Rayleigh's 1880 paper and from 
	% Michael Lugo's derivation at
	% https://gottwurfelt.com/2014/11/13/random-sums-of-sines-and-random-walks/

	nsteps = 1:100;
	for j = 1:length(nsteps)
		m(j) = mean(random_walk_demo(nsteps(j)));
	end

	% Rayleigh's formulae
	f = @(r,n) 2*r.^2./n.*exp(-r.^2./n);
	d = @(n) quad(@(r) f(r,n)/n,0,100);

	figure
	plot(nsteps,m./nsteps)
	hold on
	plot(nsteps,arrayfun(d,nsteps),'k')
	plot(nsteps,sqrt(pi*nsteps/4)./nsteps,'r--')
	legend('Numerical','Rayleigh','Lugo')

function d = random_walk_demo(nsteps,ntrials)
	if nargin < 2 || isempty(ntrials) 
		ntrials = 10000;
	end
	
	d = zeros(ntrials,1);
	for j = 1:ntrials
		d(j) = abs(sum(exp(1i*2*pi*rand(nsteps,1))));
	end