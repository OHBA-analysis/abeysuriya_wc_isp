function a = mark_zscores
	% Compute the mean and standard deviation for expected noise floor for individual connectivity
	% i.e. the values required to convert a simulation run to a z-score based on individual
	% variability in connectivity
	
	d = load('/Users/romesh/oxford_postdoc/mark_data/mark_final_connectivity_alpha.mat')

	[a.aecOrth(:,1),a.aecOrth(:,2)] = process(d.aecOrth);
	[a.plvOrth(:,1),a.plvOrth(:,2)] = process(d.plvOrth);
	[a.pli(:,1),a.pli(:,2)] = process(d.pli);

function [mc,sc] = process(m)
	% Take in N x N x S matrix (N = number of ROIs, S = number of subjects)
	% For each subject, compute the leave one out average connectivity and the
	% correlation between the individual connection matrix
	% 
	% Return - mc = mean correlation, sc = std of correlation

	flt = find(triu(ones(size(m,1)),1)); % Upper triangle entries

	c = zeros(size(m,3),1);

	for j = 1:size(m,3) % For each subject
		m2 = m;
		m2(:,:,j) = []; % Remove that subject
		m2 = mean(m2,3); % Compute the leave-one-out mean
		m3 = m(:,:,j); 
		c(j) = corr(m3(flt),m2(flt));
	end

	mc = mean(c);
	sc = std(c);