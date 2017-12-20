function format_corr(h)
	if nargin < 1 || isempty(h) 
		h = gca;
	end
	
	axis(h,'equal');
	axis(h,'tight');
	set(h,'CLim',[-1 1]);
	colormap(h,ra.plot.rwb)
	colorbar(h)