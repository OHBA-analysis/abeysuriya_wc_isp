function x = sinv_base(y,mu,sigma)
	temp = log( y./(1-y) );
	x =  mu + sigma.*temp;