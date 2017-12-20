function [c2,total_order] = reorder_matrix(c)
    % c2 = c;
    % return
    
	new_order = [29 ,14 ,8 ,33 , 34, 1  ,10 ,7 ,28 ,30 ,21 ,23 ,3  ,17 ,19 ,18 ,31 ,11 ,13 ,26 ,27 ,25 ,2  ,22 ,16 ,9  ,4, 24  ,20 ,12 ,15 ,6  ,5  ,32 ]; % Manual
	new_order = [10    12    20     4     9    24    16    28     7     1    33    30    21    23    22     3  27     2    26    17    19    18    31    25    13    11    32     5    34    29    14     8 6    15]; % GA +-1 distances
 	new_order = [    22    16    24    28     7     4    10    20    12     9     6    15     8     5    32    34 29    33    14     1    30    21    23     3    27     2    26    17    19    18    11    13 25    31]; % GA+-1 rerun
	%new_order = [ 12    10    20     4     7    28    24     9    16    22     3    23    21    30    33     1 29    14     8     6    15     5    34    32    17    19    18    11    13    25    31    26 2    27]; % GA +-2


	total_order = [new_order new_order+34]; % Reorder left and right

	% The index of total_order is the new ROI, and total_order contains the original position
	c2 = c;
	c2 = c2(total_order,:,:);
	c2 = c2(:,total_order,:);
