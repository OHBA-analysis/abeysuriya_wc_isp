function tile_cmats(cmats,varargin)
	% Plot band specific connectivity
	% cmats - cell array of n_rois x n_rois x n_subjects connectivity in each band
	% cmats will be averaged over subjects
	% 

	arg = inputParser;
	arg.addParameter('clim',[]); % Colour limits
	arg.addParameter('ncols',5); % Default number of columns
	arg.addParameter('titles',arrayfun(@(x) sprintf('Plot %d',x),1:length(cmats),'UniformOutput',false)); 
	arg.addParameter('parent',[]); 
	arg.addParameter('builtin_subplot',false);
	arg.addParameter('reorder_matrix',true); 

	arg.parse(varargin{:});

	nplots = length(cmats);
	ncols = arg.Results.ncols;
	nrows = ceil(nplots/ncols);

	if nrows == 1
		ncols = length(cmats);
	end

	if isempty(arg.Results.parent);
		parent = figure;
	else
		parent = arg.Results.parent;
	end

	for j = 1:nplots
		if arg.Results.builtin_subplot
			ax = subplot(nrows,ncols,j,'parent',parent);
		else
			ax = ra.plot.subplot(ncols,nrows,j,parent);
		end

		if arg.Results.reorder_matrix
			imagesc(diag(nan(size(ones(68),1),1))+ra.analysis.reorder_matrix(mean(cmats{j},3)),'parent',ax);
		else
			imagesc(diag(nan(size(ones(68),1),1))+mean(cmats{j},3),'parent',ax);
		end

		axis(ax,'equal')
		axis(ax,'tight')
		set(ax,'Color','k')
		%colormap(ax,'jet')
		colorbar('peer',ax)
		title(ax,arg.Results.titles{j},'interpreter','none')
		xlabel(ax,'ROI')
		ylabel(ax,'ROI')
	end


	% xlabel(ax,'ROI')
	% ylabel(ax,'ROI')
	% set(gcf,'Position',[ 26           3        1622         802])


