function x = tocol_excluding_diag(x)
	x = x(~eye(size(x)));
end