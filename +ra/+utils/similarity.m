function c = similarity(a,b)
	% Compare two correlation matrices to get their similarity
	% a and b are assumed to be correlation matrices
	% Thus the Pearson correlation is computed for the upper triangle
	% part of a and b not including the diagonal
	
	mask = logical(triu(ones(size(a)),1));
	c = corr(a(mask(:)),b(mask(:)));
