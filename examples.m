%% Wilson-Cowan simulation code
%
% This script shows how to perform simulations using the code

%% Basic usage
%
% The foundation of the code is |integrate_wc| which takes in the parameters
% for each neural population, the long-range coupling and delay matrices,
% integration duration and step size, and external input and initial
% conditions.
%
% Basic usage for the unit parameters used in the paper are shown below:
unit.mu =  [1 1];
unit.sigma =  [0.2500 0.2500];
unit.qmax =  [1 1];
unit.tau =  [0.0100 0.0200];
unit.r =  [0 0];
unit.isp_rate =  [0 0];
unit.isp_target =  [0 0];
unit.ce =  [3.5000 3.7500];
unit.ci =  [-2.5000 0];
coupl = 0;
delay = 0;
tspan = 1;
dt = 5e-5;
input(1).time = [0 tspan/2 tspan];
input(1).vals = [0.33 0.33 0.33];
input(2).time = [0 tspan/2 tspan];
input(2).vals = [0 0 0];
hist.time =  [-5.0000e-05 0];
hist.vals =  [0.4360    0.4360;...
              0.0259    0.0259];
[sol,sys] = integrate_wc( unit, coupl, delay, input, hist, tspan, dt );
figure
plot(sol.time,sol.vals(1,:));
xlabel('Time')
ylabel('Firing rate (normalized)')
title('Basic WC simulation')

%% Equivalent with WC object
% The code above is quite verbose - to simplify things, we can perform the
% same simulation using the provided Wilson-Cowan model object. First, load in
% the model parameters
u = ra.unit.deco_10

%%
% This is a struct that stores the model parameters used in the paper. Then
% construct a WilsonCowan object with 1 EI unit and set the same integration
% and sampling rate used above
wc = ra.model.WilsonCowan(u,ra.network.null(1)); 
wc.options.tspan = 1;
wc.options.tstep = 5e-5;

%%
% And run it
wc.run()
hold on
plot(wc.excitatory.time,wc.excitatory.vals,'r--')
l = legend('integrate_wc','WilsonCowan object');
set(l,'Interpreter','none')
title('Comparison of different ways to run the simulation')

%% Add noise 
% To add noise, set the noise amplitude in the unit parameters. |u.E.P| and
% |u.I.P| correspond to an external input term with mean |mu| and standard
% deviation |sigma|. So set |u.E.P.sigma| and |u.I.P.sigma| to add noise to
% the excitatory and inhibitory populations
u = ra.unit.deco_10;
u.E.P.sigma = 0.1; 
u.I.P.sigma = 0.1; 
wc = ra.model.WilsonCowan(u,ra.network.null(1)); 
wc.options.tspan = 1;
wc.options.tstep = 5e-5;
wc.run()
figure
plot(wc)
title('WC simulation with noise')

%% Changing the constant input
% If you want to change the constant input (e.g. to go from a stable fixed
% point to limit cycle oscillations), change |u.E.P.mu| and |u.I.P.mu|
u = ra.unit.deco_10;
u.E.P.mu = 0.4;
wc = ra.model.WilsonCowan(u,ra.network.null(1)); 
wc.options.tspan = 1;
wc.options.tstep = 5e-5;
wc.run()
plot(wc)
title('WC simulation with limit cycle due to extra constant input')


%% Two units
% Using |ra.network.uniform| we can construct fully connected networks with a
% given number of units. Here is one with two units/regions and no delay.
wc = ra.model.WilsonCowan(u,ra.network.uniform(2));
wc.param.coupling = [0.1 0]; % The first element is the long range global coupling, the second is E-I long range coupling (not implemented)
wc.param.velocity = inf; % Setting velocity to inf gives zero delay
wc.options.tspan = 5;
tic;wc.run();toc;
figure
plot(wc)
title('WC simulation - 2 units, no delay')

%%
% Note that the units struct is automatically expanded onto the number of
% units in the graph (if you pass in a scalar struct, if you pass in a vector
% of units they will be used directly). You can subsequently edit the
% parameters for individual regions by modifying them within the |WilsonCowan|
% object directly:
wc.units

%% Two units, with delays
% By setting the velocity to a finite value, we can include delays
wc = ra.model.WilsonCowan(u,ra.network.uniform(2));
wc.param.coupling = [0.1 0];
wc.param.velocity = 2;
wc.options.tspan = 5;
tic;wc.run();toc;
figure
plot(wc)
title('WC simulation - 2 units, with delay')

%%
% Use the |mean_delay| method to check what the mean delay is (this
% corresponds to the 'Mean Delay (s)' value shown on the plots in the paper)
wc.mean_delay


%% Cortical network
% We can simulate the cortical network by loading the provided |.mat| file.
% Note how |reorder_matrix| needs to be used to match the node ordering in the
% paper.
net = ra.network.import('data_files/cortical_network.mat');
[conns,dists] = net.netmats();
figure
imagesc(conns);
axis equal
axis tight
xlabel('ROI'); ylabel('ROI');
title('Connnectivity matrix with standard ordering used throughout code')

figure
imagesc(ra.analysis.reorder_matrix(conns));
axis equal
axis tight
xlabel('ROI'); ylabel('ROI');
title('Connnectivity matrix reordered for plotting in figures')

%%
% For consistancy, the code always uses the standard ordering, and the
% reordering step is applied at the time the figure is generated. The files
% included in the paper supplementary material are intended to be standalone
% material to accompany the paper, and thus use the ordering specified in the
% paper. This is different to the files provided in this package. If you
% are using this package, then it is important NOT to use the supplementary
% material provided with the paper, but to use the files included here that
% all use the standard ordering. All of the supplementary material in the paper
% is also provided in this package but with the standard ordering.

%% 
% The network object returned by |ra.network.import| can be passed to the
% |WilsonCowan| object
wc = ra.model.WilsonCowan(u,net);
wc.options.tspan = 1;
tic;wc.run();toc;
figure
plot(wc)
title('Cortical network simulation')

%% 
% The activity saturates quickly because global coupling of |1| is too high.
% We can reduce it to something more reasonable
wc.param.coupling = [0.04 0];

%%
% We might also want to set a particular mean delay. This can be computed
% using the distance matrix provided from |net.netmats()|, as shown below
target_delay = 10e-3; % Intended delay in seconds
mean_distance = mean(dists(logical(triu(ones(size(dists)),1))));
wc.param.velocity = mean_distance./target_delay;
wc.mean_delay % Check delay was set correctly

%%
% Now run the simulation and plot the output
tic;wc.run();toc;
figure
plot(wc)
title('Cortical network simulation, reduced coupling and shorter delay')


%% Change step size
% To change the integration step size, set |wc.options.tstep|
u = ra.unit.deco_10;
u.E.P.mu = 0.4;
wc = ra.model.WilsonCowan(u,ra.network.null(1));
wc.options.tstep = 1e-4;
wc.options.tspan = 5;
tic;wc.run();toc;
figure
plot(wc)
title('WC unit, increased simulation step size');

%% 
% Remove the initial transient - 'Burn' lets you specify a period of time to
% remove from the start and the end. Here we will remove 1s at the start to
% skip the initial transient. Note that the total simulaiton length is
% automatically adjusted so that the output still has length
% |wc.options.tspan| i.e. setting 1s of burn at the start will internally
% result in a 6s simulation, which is also reflected in the amount of time
% required to run the simulation.
tic;wc.run('burn',[1 0]);toc;
figure
plot(wc)
title('WC unit, transient removed');

%% Downsample
% Downsampling can be performed as part of running the simulation to reduce
% storage requirements. However, filtering can introduce edge artefacts at the
% start and end of the simulation. Therefore, it is useful to also specify a
% nonzero burn time at both the start and end of the simulation. In this case,
% the simulation will extend slightly beyond the requested timespan, and then
% the affected portion of the time series is discarded.
wc.options.tspan = 1;
tic;wc.run('downsample',300);toc;
figure
plot(wc)
title('Downsampled')
set(gca,'XLim',[0.9 1])

tic;wc.run('burn',[0 1],'downsample',300);toc;
hold on
plot(wc)
legend('With edge effects','With edge effects trimmed')
set(gca,'XLim',[0.9 1])

%% Set the ISP target and rate
% To enable ISP, set the ISP target and rate on the inhibitory populations in
% the unit. If you set these parameters in the excitatory population (e.g.
% |u.E.isp_target|) they will be ignored. Here, we have made tau_isp artificially
% small for the purpose of demonstrating its effect. 
u = ra.unit.deco_10;
u.I.isp_target = 0.15;
u.I.isp_rate = 5; % 1/tau_isp
u.E.P.mu = 0.35;
wc = ra.model.WilsonCowan(u,ra.network.null(1));
wc.param.coupling = [0.1 0];
wc.options.tstep = 5e-4;
wc.options.tspan = 60;
tic;wc.run('burn',[1.5 0]);toc;
figure
subplot(1,2,1)
plot(wc)
title('Neural activity')
subplot(1,2,2)
plot(wc.nsys_output.a_isp)
ylabel('c_{ie}')
title('Connection strength')

%% Rescaling firing rate
% As mentioned in the paper, the firing rate used in the model is normalized.
% This can be converted to a dimensional firing rate by multiplying the
% sigmoid by a dimensional constant, and rescaling the connection strengths
% appropriately. Changing the time constants or delays is not required, which
% means that the results can be freely rescaled with regard to firing rate
% without changing the delay axis of the results (and the connection strength
% is a single multiplicative step). Here we have an unrealistically long delay
% so that the delay effects are clearly visible (and we can thus easily verify
% that the delays do not need rescaling).
u = ra.unit.deco_10;
ics = [0.4360 0.0259 0.5497 0.4353]; % Set explicit initial conditions because these need to be rescaled as well
wc = ra.model.WilsonCowan(u,ra.network.uniform(2))
wc.param.init = ics;
wc.param.coupling = [0.14 0];
wc.param.velocity = 0.5;
wc.options.tspan = 10;
wc.run()
figure
plot(wc);
title('Normalized firing rate');

%%
% Now we select a new firing rate, and scale the connection strengths and
% initial conditions appropriately. Crucially, note that we do not change
% |u.E.tau|, |u.I.tau|, or the propagation velocity
new_qmax = 100; % Firing rate (s^{-1})
u = ra.unit.deco_10;
u.E.S.qmax = new_qmax;
u.I.S.qmax = new_qmax;
u.cee = u.cee/new_qmax;
u.cei = u.cei/new_qmax;
u.cie = u.cie/new_qmax;
wc = ra.model.WilsonCowan(u,ra.network.uniform(2))
wc.param.init = new_qmax*ics; % Change initial conditions
wc.param.coupling = [0.14/new_qmax 0]; % Change long range coupling
wc.param.velocity = 0.5; % Keep velocity the same
wc.options.tspan = 10;
wc.run()
figure
plot(wc);
ylabel('Excitatory firing rate (s^{-1})')
title('Dimensional firing rate')

%% Load a saved run and compute FC metrics
% We have included a saved simulation correspond to the optimal parameters
% shown in Fig. 5. First, load and plot the time series. Supplementary Fig.
% S5d shows a portion of the time series, so we can plot the same part of the
% signal to regenerate that figure
load(fullfile('data_files','fig_5_simulation.mat'));
figure
plot(wc)
title('Cortical network simulation');
set(gca,'XLim',[120,180],'YLim',[0 1]); % Part of the timeseries shown in Fig. S5d

%% 
% Now we show how to use the included functions to filter, orthogonalize,
% downsample, and compute connectivity metrics using AEC, PLV, and PLI. This
% regenerates the model plots shown in Fig. 5.
ts = wc.excitatory;
ts.cast('double'); % Restore double-precision
ts.trim([15 10]); % Remove start and end to ensure no edge artefacts 
[Hen,HPhase] = ts.envelope('filter',[8 13],'orthogonalize',false,'downsample',[1 NaN]);
[HenOrth,HPhaseOrth] = ts.envelope('filter',[8 13],'orthogonalize',true,'downsample',[1 NaN]);

figure
imagesc(ra.analysis.reorder_matrix(ra.analysis.aec(HenOrth))) % Note reordering step for display
colorbar
axis equal
axis tight
title('Model AEC')
figure
imagesc(ra.analysis.reorder_matrix(ra.analysis.plv(HPhaseOrth))) % Note reordering step for display
colorbar
axis equal
axis tight
title('Model PLV')
figure
imagesc(ra.analysis.reorder_matrix(ra.analysis.pli(HPhase))) % Note reordering step for display
colorbar
axis equal
axis tight
title('Model PLI')

%% 
% We can compute synchrony and metastability from the phase timecourse. The
% values below match those mentioned in the Supplementary Material (with and
% without orthogonalization).
order_ts = abs(mean(exp(1i*HPhase.vals),2));
synchrony = mean(order_ts)
metastability = std(order_ts)
order_ts = abs(mean(exp(1i*HPhaseOrth.vals),2));
synchronyOrth = mean(order_ts)
metastabilityOrth = std(order_ts)

%% 
% The subject-specific FC matrices are also provided, and can also be plotted
% (averaged over subjects). Note that the addition of |diag(nan(68,1)| sets
% the diagonal of the matrix to |NaN| so that the colour range automatically
% ignores the self connections.
d = load(fullfile('data_files','data_fc.mat'));

figure
imagesc(diag(nan(68,1))+ra.analysis.reorder_matrix(mean(d.aecOrth,3))) % Note reordering step for display
colorbar
axis equal
axis tight
title('Data AEC')

figure
imagesc(diag(nan(68,1))+ra.analysis.reorder_matrix(mean(d.plvOrth,3))) % Note reordering step for display
colorbar
axis equal
axis tight
title('Data PLV')

figure
imagesc(diag(nan(68,1))+ra.analysis.reorder_matrix(mean(d.pli,3))) % Note reordering step for display
colorbar
axis equal
axis tight
title('Data PLI')

%% Reproducing Fig. 5 from scratch
% If we wanted to reproduce the simulations in Fig. 5 from scratch, we could
% construct it as follows. See |sweeps/deco_rk4_isp.m| and
% |+ra/+sweep/simulate_plasticity.m|
%
%  u = ra.unit.deco_10;
%  noise_mu = [0.31 0];
%  noise_sigma = [1e-2 1e-2];
%  u.E.P = struct('mu',noise_mu(1),'sigma',noise_sigma(1),'name','gaussian');
%  u.I.P = struct('mu',noise_mu(2),'sigma',noise_sigma(2),'name','gaussian');
%  u.I.isp_target = 0.15;
%  
%  net = ra.network.import('stam_cortical');
%  [c,d] = net.netmats;
%  mean_distance = mean(d(logical(triu(ones(size(d)),1))));
%  min_distance = min(d(logical(triu(ones(size(d)),1))));
%  mean_delays = 1e-3*(9);
%  velocity = mean_distance./mean_delays;
%  
%  wc = ra.model.WilsonCowan(u,net);
%  wc.param.coupling = [0.14 0];
%  wc.param.velocity = velocity;
%  wc.options.tstep = 1e-4;  
%  wc.mx_helper = 'wc_network_rk4';
%  wc.param.round_delays = false;
%  [wc,isp_ts,activity_ts] = ra.wc.ramped_plasticity(wc,50,30) % Long ISP run (1500s)
%  wc = ra.wc.long_run(wc,25,21,300); % Output run for analyzing activity
%
% This would be expected to take on the order of hours to simulate.
