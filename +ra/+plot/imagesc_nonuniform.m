function pc = imagesc_nonuniform(x,y,data)
	% Like imagesc, but supports nonuniform X and Y
	% First, pad the data
	data(end+1,:) = NaN;
	data(:,end+1) = NaN;

	% Now we assume the first row of data corresponds to the first value of y
	% and first column of data corresponds to the first value of x
	% same as what comes out of meshgrid (NOT ndgrid)
	
	if x(end) < x(1)
		x = x(end:-1:1);
		data = data(:,end:-1:1);
	end

	if y(end) < y(1)
		y = y(end:-1:1);
		data = data(end:-1:1,:);
	end
	
	x(end+1) = x(end)+(x(end)-x(end-1))/2;
	y(end+1) = y(end)+(y(end)-y(end-1))/2;
	pc = pcolor(x,y,data);
	set(pc,'LineStyle','none')

	% Now the display is correct, but the axis labels and ticks need to be adjusted
	set(gca,'XTick',x(1:end-1)+diff(x)/2,'XTickLabel',arrayfun(@(x) sprintf('%.2g',x),x(1:end-1),'UniformOutput',false))
	set(gca,'YTick',y(1:end-1)+diff(y)/2,'YTickLabel',arrayfun(@(x) sprintf('%.2g',x),y(1:end-1),'UniformOutput',false))

