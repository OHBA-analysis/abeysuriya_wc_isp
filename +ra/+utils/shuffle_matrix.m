function [y,idx] = shuffle_matrix(x,idx)
	if nargin < 2 || isempty(idx) 
		idx = find(triu(logical(ones(size(x))),1));
	end

	order = randperm(length(idx));
	y = zeros(size(x));
	y(idx(order)) = x(idx);
	y = y + y.';