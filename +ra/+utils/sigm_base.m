function y = sigm_base(x,mu,sigma)
	x = (x-mu) ./ sigma ; % rescale
	x = exp(-x);
	y =  1./(1+x);
