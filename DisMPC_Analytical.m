%% Distributed MPC Algorithm for Spacecraft Rendezvous
% =========================================================================
% AA277  |  Luke Neise, Samuel Low, Michael Ying, Tamas Kis

clc; clear all; close all;
dt = 60;           % Dynamics time step
a = 6725000;       % Semi-major axis [m]
N = 5;            % MPC prediction horizon size
u_lb = -3.0;       % ΔV lower bound [m/s]
u_ub =  3.0;       % ΔV upper bound [m/s]
Q = 0.1*eye(6);    % State error cost matrix
R = 0.001*eye(3);  % Control effort cost matrix
P = N^2*eye(6);    % Terminal state error cost matrix

% Because the semi-major axis is a much larger numerical quantity than all
% the other angular variables, the optimization weighting should scale it
% down by the inverse arc length constant for relative weighting.
Q(1,1) = 0.2*N / (0.0175 * a); 
P(1,1) = 0.2*N / (0.0175 * a);



%% Exterior Online Simulation
% =========================================================================
% Set up our spacecraft and reference
% Elements: [ a, ex, ey, inc, argp, nu ]
xk1 = [a+25000,  0.002,  0.002, deg2rad(45.05), deg2rad(45.05), 1*pi/180]'; 
xk2 = [a-25000, -0.002, -0.002, deg2rad(45.00), deg2rad(45.00),-1*pi/180]';

% Setup some arrays for plotting the final results.
xF1a = [];
xF2a = [];
uF1a = [];
uF2a = [];
dv1a = [0];
dv2a = [0];

% Initialize some simulation time
duration = 86400;
timeaxis = linspace( 0, duration, (duration/dt)+1 );

% Initialize arrays for run time
t_runtime = zeros( 1, length(timeaxis)-1 );
t_runtimesum = zeros( 1, length(timeaxis)-1 );

% Main for loop.
for k = 0 : 1 : round(duration/dt)

    % Run the main DMPC program
    tic;
    [xF1, xF2, uF1, uF2] = run_dist_MPC( xk1, xk2, dt, ...
                                         N, Q, R, P, ...
                                         u_lb, u_ub );

    % Update the states in external simulation
    xF1a(:,k+1) = xF1(:,1); % Take only the first element.
    xF2a(:,k+1) = xF2(:,1); % Take only the first element.
    uF1a(:,k+1) = uF1(:,1); % Take only the first element.
    uF2a(:,k+1) = uF2(:,1); % Take only the first element.
    xk1 = nonlinear_dynamics( xF1(:,1), uF1(:,1), dt );
    xk2 = nonlinear_dynamics( xF2(:,1), uF2(:,1), dt );
    
    % Update the Delta V values for SC1 and SC2
    dv1 = sum(abs(uF1(:,1)));
    dv1a(end+1) = dv1a(end) + dv1;
    dv2 = sum(abs(uF2(:,1)));
    dv2a(end+1) = dv2a(end) + dv2;

    % Update the run-time
    t_runtime(k+1) = toc;
    if k+1 > 1
        t_runtimesum(k+1) = t_runtimesum(k) + t_runtime(k+1);
    else
        t_runtimesum(k+1) = t_runtime(k+1);
    end
end

% For distributed MPC, backpropagate the states. This step is not needed
% for decentralized MPC, where you can simply refer to the states of the 
% reference from before.
mu = 3.986004415e14;
xFR = 0.5 * ( xF1a(:,end) + xF2a(:,end) );
xFRa = [];
xFRa(:,1) = xFR;
nR = -1*sqrt(mu/(xFR(1))^3);
for k = 1 : 1 : round(duration/dt)
    xFRa(:,k+1) = xFRa(:,k);
    xFRa(6,k+1) = xFRa(6,k) + (nR*dt);
end

% Flip the reference orbit because we backpropagated it
xFRa = flip(xFRa,2); 

% Loop through all the states again now that we have the final orbit and
% perform the conversion to RTN states
RTN_pos1 = [];
RTN_pos2 = [];
RTN_vel1 = [];
RTN_vel2 = [];
for k = 1 : 1 : round(duration/dt)
    oe1 = xF1a(:,k);
    oe2 = xF2a(:,k);
    oeR = xFRa(:,k);
    % RTN1 = elements_to_RTN( oeR, oe1 );
    % RTN2 = elements_to_RTN( oeR, oe2 );
    RTN1 = elements_to_RTN( oe1, oe2 );
    RTN2 = RTN1;
    RTN_pos1(:,k) = RTN1(1:3)';
    RTN_pos2(:,k) = RTN2(1:3)';
    RTN_vel1(:,k) = RTN1(4:6)';
    RTN_vel2(:,k) = RTN2(4:6)';
end

% converts time to hours
timeaxis = timeaxis/3600;

% converts semi-major axes to kilometers
xF1a(1,:) = xF1a(1,:)/1000;
xF2a(1,:) = xF2a(1,:)/1000;

% converts angles to degrees
xF1a(4,:) = xF1a(4,:)*(180/pi);
xF2a(4,:) = xF2a(4,:)*(180/pi);
xF1a(5,:) = xF1a(5,:)*(180/pi);
xF2a(5,:) = xF2a(5,:)*(180/pi);
xF1a(6,:) = xF1a(6,:)*(180/pi);
xF2a(6,:) = xF2a(6,:)*(180/pi);

% converts RTN positions to kilometers
RTN_pos1 = RTN_pos1/1000;


%% Saving data

% time [h]
t = timeaxis;

% semi-major axes [km]
a1 = xF1a(1,:);
a2 = xF2a(1,:);

% x-component of eccentricity vector [-]
ex1 = xF1a(2,:);
ex2 = xF2a(2,:);

% y-component of eccentricity vector [-]
ey1 = xF1a(3,:);
ey2 = xF2a(3,:);

% inclination [deg]
i1 = xF1a(4,:);
i2 = xF2a(4,:);

% RAAN [deg]
Om1 = xF1a(5,:);
Om2 = xF2a(5,:);

% argument of latitude [deg]
u1 = xF1a(6,:);
u2 = xF2a(6,:);

% control effort [m/s]
ctrl1 = uF1a;
ctrl2 = uF2a;

% cumulative Delta-V [m/s]
dV1 = dv1a;
dV2 = dv2a;

% RTN positions [km]
RTN_2wrt1 = RTN_pos1;

% trims all arrays so they are same length (not sure why they are
% different?) and packages them into a structure
DMPC_analytical.t = t(1:1440);
DMPC_analytical.a1 = a1(1:1440);
DMPC_analytical.a2 = a2(1:1440);
DMPC_analytical.ex1 = ex1(1:1440);
DMPC_analytical.ex2 = ex2(1:1440);
DMPC_analytical.ey1 = ey1(1:1440);
DMPC_analytical.ey2 = ey2(1:1440);
DMPC_analytical.i1 = i1(1:1440);
DMPC_analytical.i2 = i2(1:1440);
DMPC_analytical.Om1 = Om1(1:1440);
DMPC_analytical.Om2 = Om2(1:1440);
DMPC_analytical.u1 = u1(1:1440);
DMPC_analytical.u2 = u2(1:1440);
DMPC_analytical.ctrl1 = ctrl1(:,1:1440);
DMPC_analytical.ctrl2 = ctrl2(:,1:1440);
DMPC_analytical.dV1 = dV1(1:1440);
DMPC_analytical.dV2 = dV2(1:1440);
DMPC_analytical.RTN_2wrt1 = RTN_2wrt1(:,1:1440);

% saves simulation data
save('data/DMPC_analytical.mat','DMPC_analytical');



%% Primary Distributed MPC Optimization Process (Analytical).
% =========================================================================

function [xf1, xf2, uf1, uf2] = run_dist_MPC( xk1, xk2, dt, ...
                                              N, Q, R, P, u_lb, u_ub )
    % Input: xk1, xk2, are 6xN matrices detailing the states of SC1/2
    % Initialize MPC global and local variables.
    u10 = zeros(3,N); % Initialize RTN control ΔV for SC1 [m/s]
    u20 = zeros(3,N); % Initialize RTN control ΔV for SC2 [m/s]
    uf1 = u10;
    uf2 = u20;
    xf1 = xk1(:,1); 
    xf2 = xk2(:,1);
    mu = 3.986004415e14;
    Rinv = inv(R);

    % Re-initialize (re-update) the initial state to length 1.
    xf1 = xf1(:,1); % Take just the first vector.
    xf2 = xf2(:,1); % Take just the first vector.

    % Populate the states of both spacecraft first.
    for i = 1:N
        xf1(:,i+1) = nonlinear_dynamics( xf1(:,i), uf1(:,i), dt );
        xf2(:,i+1) = nonlinear_dynamics( xf2(:,i), uf2(:,i), dt );
    end

    % (Re)-populate states based on control from previous iteration.
    for i = 1:N
        
        % This is supposed to be a while loop to check for tolerance but
        % got lazy and decided to just 5 iterations for iterative DMPC. 
        % The variable uf1 at the end of the loop will be fed into the
        % first two lines xf1(:,i+1) and xf2(:,i+1) (the iterative part)
        for q = 1:5

            xf1(:,i+1) = nonlinear_dynamics( xf1(:,i), uf1(:,i), dt );
            xf2(:,i+1) = nonlinear_dynamics( xf2(:,i), uf2(:,i), dt );
            
            % Now compute the B matrix in the state transition for SC1.
            sma1 = xf1(1,i);
            inc1 = xf1(4,i);
            lat1 = xf1(6,i);
            n1 = sqrt(mu/(sma1)^3);
            B1 = [ 0   n1   0 ;
                   sin(lat1)/(sma1*n1)   2*cos(lat1)/(sma1*n1)   0 ;
                   -cos(lat1)/(sma1*n1)   2*sin(lat1)/(sma1*n1)   0 ;
                   0   0   cos(lat1)/(sma1*n1) ;
                   0   0   sin(lat1)/(sma1*n1*sin(inc1)) ;
                   -2/(sma1*n1)   0   -cot(inc1)*sin(lat1)/(sma1*n1) ];
    
            % Now compute the B matrix in the state transition for SC2.
            sma2 = xf2(1,i);
            inc2 = xf2(4,i);
            lat2 = xf2(6,i);
            n2 = sqrt(mu/(sma2)^3);
            B2 = [ 0   n2   0 ;
                   sin(lat2)/(sma2*n2)   2*cos(lat2)/(sma2*n2)   0 ;
                   -cos(lat2)/(sma2*n2)   2*sin(lat2)/(sma2*n2)   0 ;
                   0   0   cos(lat2)/(sma2*n2) ;
                   0   0   sin(lat2)/(sma2*n2*sin(inc2)) ;
                   -2/(sma2*n2)   0   -cot(inc2)*sin(lat2)/(sma2*n2) ];

            % Compute the horizon-stacked transitions for SC1 and SC2.
            A1 = eye(6);
            A2 = eye(6);
            A1(6,1) = (-3*n1*dt/(2*sma1));
            A2(6,1) = (-3*n2*dt/(2*sma2));
            A1stack = A1;
            A2stack = A2;
            for m = 1:(N-i)
                A1stack = A1stack * A1;
                A2stack = A2stack * A2;
            end
            
            % Compute the optimal control vector for SC1
            x1_0_term = -1 * B1' * Q * ( xf1(:,i+1) - xf2(:,i+1) );
            x1_F_term = -1 * B1' * A1stack * P * (xf1(:,end)-xf2(:,end));
            uf1(:,i) = x1_0_term + x1_F_term;
            uf1(:,i) = Rinv * uf1(:,i);
            
            % Compute the optimal control vector for SC2
            x2_0_term = -1 * B2' * Q * ( xf2(:,i+1) - xf1(:,i+1) );
            x2_F_term = -1 * B2' * A2stack * P * (xf2(:,end)-xf1(:,end));
            uf2(:,i) = x2_0_term + x2_F_term;
            uf2(:,i) = Rinv * uf2(:,i);

            % Take the saturation if the control input exceeds.
            uf1 = max( uf1, u_lb );
            uf1 = min( uf1, u_ub );
            uf2 = max( uf2, u_lb );
            uf2 = min( uf2, u_ub );
            
            % Recompute next state
            xf1(:,i+1) = nonlinear_dynamics( xf1(:,i), uf1(:,i), dt );
            xf2(:,i+1) = nonlinear_dynamics( xf2(:,i), uf2(:,i), dt );
        
        end
    end
end