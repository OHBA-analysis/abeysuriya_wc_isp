function param = deco_10(scale)
    % Sub-threshold Deco unit
    % Oscillations at around 11Hz
    % PARAMETERS FROM DECO - see below
    %
    % Deco uses sigmoid

    % 1/sigma = 0.1 -> sigma = 10
    % mu = 40
    % qmax = 100
    % tau_e = 1
    % tau_i = 0.2
    % P_e = 10
    % c_ee = 1.4
    % c_ei = 1.5
    % c_ii = 0
    % c_ie = 1
    % First, rescale the sigmoid to unit activation, we just divide qmax by 100 and leave everything
    % else unchanged - and multiply the coupling strengths by 100 (don't forget x value in sigmoid is not in the same units as firing)
    % Next, we rescale the sigmoid mean. Scaling just means multiplying both by a constant
    % So we divide both by 40 to get 

    % sigma = 0.25
    % mu = 1

    % And in doing so, we get 

    % c_ee = 1.4/4 = 3.50
    % c_ei = 1.5/40 = 3.75
    % c_ii = 0
    % c_ie = 1/40 = 2.5

    % Lastly, we set the time constants. Note that we need to swap tau and tau' compared to the paper in order to have a hopf bifurcation. We seee

    % tau_i = 1/100 -> 0.01
    % tau_e = 0.2/100 -> 0.002
    % This has a bifurcation at 28.51 Hz
    % u.I.tau = 0.006
    %This gives a bifurcation at 40 Hz

    if nargin < 1 || isempty(scale) 
        scale = 1;
    end
   
    % Excitatory subpop
    E.rp    = 0; 
    E.tau   = 0.01; %2.75e-3; % Based on getting the LLE 

    E.S.name = 'logistic';
    E.S.mu    = 1*scale;
    E.S.sigma = 0.25*scale;

    E.P.name = 'gaussian';
    E.P.mu = 0.33;
    E.P.sigma = 0;

    % Inhibitory subpop
    I.rp    = 0; 
    I.tau   = 2*E.tau; 

    I.S.name = 'logistic';
    I.S.mu    = 1*scale;
    I.S.sigma = 0.25*scale;

    I.P.name = 'gaussian';
    I.P.mu = 0;
    I.P.sigma = 0;


    % local couplings
    param.cee =  3.5*scale;
    param.cie =  -2.5*scale;
    param.cei =   3.75*scale;
    param.cii =   0;
    
    E.isp_rate = 0;
    E.isp_target = 0;
    E.isp_stop_time = 0;

    E.hip_rate = 0;
    E.hip_target = 0;
    E.hip_stop_time = 0;

    I.isp_rate = 0;
    I.isp_target = 0;
    I.isp_stop_time = 0;

    I.hip_rate = 0;
    I.hip_target = 0;
    I.hip_stop_time = 0;


    param.E = E;
    param.I = I;

