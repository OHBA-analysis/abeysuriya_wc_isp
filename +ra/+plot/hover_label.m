function hover_label(p)
	% Take in parcellation, change title on hover
	ax = gca;
	assert(diff(get(ax,'XLim')) == p.n_parcels )
	assert(diff(get(ax,'YLim')) == p.n_parcels )


	[~,total_order] = ra.analysis.reorder_matrix(ones(p.n_parcels));
	labels = p.labels(total_order);

	set(gcf,'Units','pixels');
	set(ax,'Units','normalized');

	set(gcf,'WindowButtonMotionFcn',@(a,b) changetitle(a,ax,labels))

	set(gca,'XTick',[1:2:p.n_parcels],'XTickLabel',labels(1:2:end),'XTickLabelRotation',90)

	set(gca,'YTick',[2:2:p.n_parcels],'YTickLabel',labels(2:2:end));

function within = is_within(C,pos)
	if C(1) < pos(1) || C(1) > pos(1)+pos(3) || C(2) < pos(2) || C(2) > pos(2)+pos(4)
		within = false;
	else
		within = true;
	end

function changetitle(a,ax,labels)

	if ~is_within(get(a,'CurrentPoint'),getpixelposition(ax))
		title(ax,'');
		return
	else
		% Get parcel index
		C = get(ax,'CurrentPoint');
		C = round(C(1,1:2));
		l = labels(C);
		title(sprintf('%s -> %s',l{1},l{2}))
	end

