function h = vector_imagesc(x,y,data,varargin)
	% Like imagesc, but supports nonuniform X and Y
	% First, pad the data
	if ~isvector(x)
		data = x;
		x = 1:size(data,1);
		y = 1:size(data,2);
	end
	
	data(end+1,:) = NaN;
	data(:,end+1) = NaN;

	% Flip the dimensions in cases where the X and Y arrays are descending instead of ascending
	if x(end) < x(1)
		x = x(end:-1:1);
		data = data(:,end:-1:1);
	end

	if y(end) < y(1)
		y = y(end:-1:1);
		data = data(end:-1:1,:);
	end
	
	% Add a dummy point at the end
   
	x(end+1) = x(end)+(x(end)-x(end-1));
	y(end+1) = y(end)+(y(end)-y(end-1));

	% Move the whole thing half a space
	x = x - (x(2)-x(1))/2;
	y = y - (y(2)-y(1))/2;

	cl = [min(data(:)),max(data(:))];
	data(isnan(data)) = 0;

	h = pcolor(x,y,data,varargin{:});
	set(h,'LineStyle','none')
	ax = get(h,'Parent');
	set(ax,'XLim',[min(x) max(x)]);
	set(ax,'yLim',[min(y) max(y)]);
	set(ax,'YDir','reverse','Layer','top')
	set(ax,'CLim',cl);
