function pop_out()
	ax = gca;
	f2=figure;
	h2 = copyobj(ax,f2)

	set(h2,'Position',[0.1300 0.1100 0.7750 0.8150])

	if ~isempty(findobj(f2,'Type','image'))
		colorbar
	end

	set(findobj(gcf),'HitTest','on');
	
