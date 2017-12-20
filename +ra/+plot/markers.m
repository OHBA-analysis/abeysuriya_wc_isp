function markers(bg,ax)
	% Set axes callback to put an incremented marker on the plot for each click
	% e.g. plot(rand(5,1)); ra.draw_markers; <click on the plot a couple of times>
	%
	% This will switch the plot into marker mode i.e.
	% hittest for all children is turned off, so you can 
	% click anywhere
	if nargin < 2 || isempty(ax) 
		ax = gca;
	end
	
	if nargin < 1 || isempty(bg) 
		bg = 'none';
	end
	
	try
		c = findobj(ax,'HitTest','on');
		set(c,'HitTest','off');
	end
	set(ax,'HitTest','on')

	set(ax,'UserData',double('a'));
	set(ax,'ButtonDownFcn',@(a,b,c) add_text(a,b,ax,bg))


function add_text(a,b,ax,bg)
	clicked = b.IntersectionPoint(1:2);
	hold(ax,'on')
	text(clicked(1),clicked(2),sprintf('(%s)',get(ax,'UserData')),'Color','r','Fontsize',13,'HorizontalAlignment','center','VerticalAlignment','middle','FontWeight','bold','background',bg);
	set(ax,'UserData',get(ax,'UserData')+1);
