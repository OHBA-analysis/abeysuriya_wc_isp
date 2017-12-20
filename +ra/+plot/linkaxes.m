function linkaxes(x,y)
	if nargin < 1 || isempty(x) 
		x = true;
	end

	if nargin < 2 || isempty(y) 
		y = true;
	end
	
	
	figs = get(0,'Children');
	ax = [];
	for j = 1:length(figs)
		ax = [ax findobj(figs(j),'Type','axes')];
	end

	if x && y
		linkobj = linkprop(ax,{'XLim','YLim','View'})
	elseif x
		linkobj = linkprop(ax,{'XLim','View'})
	elseif y
		linkobj = linkprop(ax,{'YLim','View'})
	else
		linkobj = linkprop(ax,{'View'})
	end

	set(figs(1),'UserData',linkobj)