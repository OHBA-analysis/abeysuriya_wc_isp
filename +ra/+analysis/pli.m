function pli_out = pli(HPhase)
	% Take in phase timecourses and compute phase lag index
	%
	% Example usage
	%
	% 	[~,HPhase] = ts.envelope('filter',[low high])
	% 	ra.analysis.pli(HPhase)
	%
	% Artifacts should have been removed first e.g.
	%
	% 	ra.analysis.pli(HPhase.select_times(clean))

	assert(all(isfinite(HPhase.vals(:)))); % Any NaNs should have been removed

	pli_out = zeros(HPhase.n_signals,HPhase.n_signals,length(HPhase));

	for j = 1:length(HPhase)
		[a,b] = meshgrid(1:HPhase.n_signals);
		a = triu(a,1);
		b = triu(b,1);
		a = a(:);
		b = b(:);
		idx = find(a~=0 & b~=0);

		c = zeros(HPhase.n_signals,HPhase.n_signals);
		HPhase_diffs = zeros(HPhase.n_times,length(idx));

		for l = 1:length(idx)
			HPhase_diffs(:,l) = HPhase.vals(:,a(idx(l)))-HPhase.vals(:,b(idx(l)));
		end

		tmp = abs(mean(sign(sin(HPhase_diffs)),1));

		for l = 1:length(idx)
			c(a(idx(l)),b(idx(l))) = tmp(l);
		end

		c = tril(c)+tril(c)';

		pli_out(:,:,j) = c;

	end

	