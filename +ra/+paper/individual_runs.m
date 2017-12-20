function out = individual_runs

	out.isp = fullfile(startup.get_rootdir,'romesh_nsys','sweeps','deco_rk4_isp','results_part_262.mat'); % 0.1400, mean_delay = 0.0080
	out.noisp = fullfile(startup.get_rootdir,'romesh_nsys','sweeps','deco_rk4','results_part_76.mat'); 

	%out.noisp = load(fullfile(startup.get_rootdir,'romesh_nsys','sweeps','deco_rk4','results_part_97.mat')); % most resembles ISP case but PLV and PLI aren't as good.
	%out.noisp = load(fullfile(startup.get_rootdir,'romesh_nsys','sweeps','deco_rk4','results_part_53.mat')); % Part 53 has good AEC but poor PLI and very different timeseries

	% Extra runs
	out.preisp = fullfile(startup.get_rootdir,'romesh_nsys','sweeps','deco_rk4','results_part_261.mat'); % Should match out.isp run number
	
	% These are for use in isp_convergence figure, where out.isp is the mid synchrony case
	out.high_synchrony = fullfile(startup.get_rootdir,'romesh_nsys','sweeps','deco_rk4_isp','results_part_444.mat');   % 0.2300 , delay  0.0020
	out.low_synchrony = fullfile(startup.get_rootdir,'romesh_nsys','sweeps','deco_rk4_isp','results_part_78.mat'); % 0.05, 0.014

	% These are for use in the isp_
	out.high_isp_matching = fullfile(startup.get_rootdir,'romesh_nsys','sweeps','deco_rk4_isp_target_30_hd','results_part_262.mat'); % 0.05, 0.014
	out.high_isp_worst = fullfile(startup.get_rootdir,'romesh_nsys','sweeps','deco_rk4_isp_target_30_hd','results_part_103.mat'); % 0.05, 0.014

	out.noise_0x = fullfile(startup.get_rootdir,'romesh_nsys','sweeps','rev2_noise','deco_rk4_isp_noise_test','results_part_1.mat');
	out.noise_5x = fullfile(startup.get_rootdir,'romesh_nsys','sweeps','new','deco_rk4_isp_noise_test_20170913_1','results_part_1.mat'); 
	out.noise_10x = fullfile(startup.get_rootdir,'romesh_nsys','sweeps','new','deco_rk4_isp_noise_test_20170913_2','results_part_1.mat'); 
	out.noise_20x = fullfile(startup.get_rootdir,'romesh_nsys','sweeps','new','deco_rk4_isp_noise_test_20170913_3','results_part_1.mat'); 
