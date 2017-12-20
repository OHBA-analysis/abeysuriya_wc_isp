function mapping = map_times(t1,t2)
	% Given a time vector like
	% t1 = [1 2 3 4 5 6 7 8 9 10]
	% And downsampled times like
	% t2 = [1 4 8]
	% Find which index in t2 is closest to t1
	% e.g
	% mapping = [1 1 2 2 2 3 3 3 3 3]

	assert(isvector(t1) && isvector(t2),'Inputs must be vectors');
	assert(issorted(t1) && issorted(t2),'Inputs must be sorted');

	ptr = 1;
	mapping = zeros(size(t1));

	for j = 1:length(t1)
		d = abs(t1(j) - t2([ptr ptr+1]));
		if d(1) >= d(2)
			ptr = ptr + 1;
		end

		if ptr == length(t2)
			mapping(j:end) = ptr;
			break
		else
			mapping(j) = ptr;
		end

	end
