function ax = subplot(nx,ny,idx,parent)
	if nargin < 4 || isempty(parent) 
		parent = gcf;
	end
	
	x_margin = 0.02;
	y_margin = 0.05;
	x_pad = 0.01;
	y_pad = 0.06;

	x_width = (1-2*x_margin-(nx-1)*x_pad)/nx;
	y_width = (1-2*y_margin-(ny-1)*y_pad)/ny;

	[a,b] = ind2sub([nx ny],idx);
	b = ny-b+1;

	x_pos = x_margin+(x_width+x_pad)*(a-1);
	y_pos = y_margin+(y_width+y_pad)*(b-1);
	ax = axes('Position',[x_pos y_pos x_width y_width],'Units','normalized','parent',parent);
end