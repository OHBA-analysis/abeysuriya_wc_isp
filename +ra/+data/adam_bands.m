function [band_names,band_frequencies] = adam_bands(selector)
	if nargin < 1 || isempty(selector) 
		selector = 1;
	end
	
	switch selector
		case 1
			band_names = {'c_delta','theta','c_low_alpha','alpha','c_wide_beta','beta','c_beta_gamma','low_gamma','mid_gamma','high_gamma','bb_1_85','bb_1_30','bb_4_30'};
			band_frequencies = {[2 6] ,[4 8] ,[6 10.5] ,[8 13] ,[10.5 21.5] ,[13 30] ,[21.5 39] ,[30 48] ,[39 66] ,[52 80],[1 85],[1 30],[4 30]};
		case 2
			band_names =  {'2_6','4_8','6_10','8_12','10_14','12_16','14_18','16_20','18_22','20_24','22_26','24_28','26_30'};
			band_frequencies = {[2 6] ,[4 8] ,[6 10] ,[8 12] ,[10 14] ,[12 16] ,[14 18] ,[16 20] ,[18 22] ,[20 24] ,[22 26] ,[24 28] ,[26 30]};
		case 3
			centers = 4:2:80;
			for j = 1:length(centers)
				band_names{j} = sprintf('%d_%d',centers(j)-2,centers(j)+2);
				band_frequencies{j} = [-2 2]+centers(j);
			end
		case 4
			% Let's go in wider windows - instead of 4Hz, consider 5Hz 
			error('Not yet')
			band_names =  {'2_6','4_8','6_10','8_12','10_14','12_16','14_18','16_20','18_22','20_24','22_26','24_28','26_30'};
			band_frequencies = {[2 6] ,[4 8] ,[6 10] ,[8 12] ,[10 14] ,[12 16] ,[14 18] ,[16 20] ,[18 22] ,[20 24] ,[22 26] ,[24 28] ,[26 30]};
		otherwise
			band_names = {'c_delta','theta','c_low_alpha','alpha','c_wide_beta','beta','c_beta_gamma','bb_1_30','bb_4_30'};
			band_frequencies = {[2 6] ,[4 8] ,[6 10.5] ,[8 13] ,[10.5 21.5] ,[13 30] ,[21.5 39],[1 30],[4 30]};
		end
	end

	













