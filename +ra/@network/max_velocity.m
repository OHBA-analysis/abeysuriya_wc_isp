function v = max_velocity(self,dt)
	% Given a timestep, determine the largest velocity compatible with a certain minimum
	% delay. Specifically, the velocity should NEVER be large enough that it is faster
	% than the smallest timestep. 

	if nargin < 1 || isempty(dt) 
		error('You must specify the timestep to calculate the max permitted velocity')
	end

	min_dist = min(self.edges.Distance(self.edges.Distance>0)); % Smallest nonzero distance
	
	for j = 1:10
	 	min_delay = dt*j;
	 	max_velocity = min_dist/min_delay;
	 	fprintf('Min # steps = %d, Max velocity = %.2f\n',j,max_velocity);
	end

	v = floor(min_dist/(dt*2)); % Return velocity based on threshold of 2 steps as checked by ra.model.WilsonCowan