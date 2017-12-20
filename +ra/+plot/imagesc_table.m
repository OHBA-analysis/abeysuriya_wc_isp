function imagesc_table(data,clim)
	% Plot an imagesc map and superimpose text labels with the values on the plot
	
	if nargin < 2 || isempty(clim) 
		clim = [];
	end
	
	figure
	imagesc(data)
	axis equal
	axis tight
	if ~isempty(clim)
		set(gca,'CLim',clim);
	end
	colorbar

	textvec = []
	cl = get(gca,'CLim');
	pc = @(x) (x-cl(1))./(diff(cl)); % 
	for j = 1:size(data,1)
		for k = 1:size(data,2)
			if pc(data(j,k)) < 0.75
				color = 'w';
			else
				color = 'k';
			end
			textvec(end+1) =  text(j,k,sprintf('%.3f',data(j,k)),'Color',color,'VerticalAlignment','middle','HorizontalAlignment','center','Fontsize',16)
		end
	end
