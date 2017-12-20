function plv_out = plv(HPhase)
	% Take in phase timecourses and compute phase locking value
	%
	% Example usage
	%
	% 	[~,HPhase] = ts.envelope('filter',[low high])
	% 	ra.analysis.plv(HPhase)
	%
	% Artifacts should have been removed first e.g.
	%
	% 	ra.analysis.plv(HPhase.select_times(clean))
	
	assert(all(isfinite(HPhase.vals(:)))); % Any NaNs should have been removed

	plv_out = zeros(HPhase.n_signals,HPhase.n_signals,length(HPhase));

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

		tmp = abs(mean(exp(1i*HPhase_diffs),1));

		for l = 1:length(idx)
			c(a(idx(l)),b(idx(l))) = tmp(l);
		end

		c = tril(c)+tril(c)';

		plv_out(:,:,j) = c;

	end

	