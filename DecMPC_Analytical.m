%% Decentralized MPC Algorithm for Spacecraft Rendezvous
% =========================================================================
% AA277  |  Luke Neise, Samuel Low, Michael Ying, Tamas Kis

clc; clear all; close all;
dt = 60;           % Dynamics time step
a = 6700000;       % Semi-major axis [m]
N = 50;            % MPC prediction horizon size
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
xk1 = [a+25000,  0.002,  0.002, deg2rad(90.05), deg2rad(90.05), 0]'; 
xk2 = [a-25000, -0.002, -0.002, deg2rad(89.95), deg2rad(89.95), 0]';
xkR = [a      , 0, 0, deg2rad(90.00), deg2rad(90.00), 0]';

% Setup some arrays for plotting the final results.
xF1a = [];
xF2a = [];
xFRa = [];
uF1a = [];
uF2a = [];
dv1a = [0];
dv2a = [0];

% Initialize some simulation time
duration = 12000;
timeaxis = linspace( 0, duration, (duration/dt)+1 );

for k = 0 : 1 : round(duration/dt)
    [xF1, xF2, xFR, uF1, uF2] = run_DecMPC( xk1, xk2, xkR, dt, ...
                                            N, Q, R, P, ...
                                            u_lb, u_ub );
    xF1a(:,k+1) = xF1(:,1); % Take only the first element.
    xF2a(:,k+1) = xF2(:,1); % Take only the first element.
    xFRa(:,k+1) = xFR(:,1); % Take only the first element.
    uF1a(:,k+1) = uF1(:,1); % Take only the first element.
    uF2a(:,k+1) = uF2(:,1); % Take only the first element.
    xk1 = nonlinear_dynamics( xF1(:,1), uF1(:,1), dt );
    xk2 = nonlinear_dynamics( xF2(:,1), uF2(:,1), dt );
    xkR = nonlinear_dynamics( xFR(:,1), zeros(3,1), dt );
    dv1 = sum(abs(uF1(:,1)));
    dv1a(end+1) = dv1a(end) + dv1;
    dv2 = sum(abs(uF2(:,1)));
    dv2a(end+1) = dv2a(end) + dv2;
end



%% Plotting portion of the code.
% =========================================================================

figure(1) % Semi-major axis (m)
plot( timeaxis, xF1a(1,:), LineWidth=1.25 )
hold('on')
grid('on')
plot( timeaxis, xF2a(1,:), LineWidth=1.25 )
plot( timeaxis, xFRa(1,:), 'k--', LineWidth=1.25 )
title('Semi-major axis of SC1 and SC2')
legend('SC1','SC2','Reference')

figure(2) % Eccentricity vector X (unitless)
plot( timeaxis, xF1a(2,:), LineWidth=1.25 )
hold('on')
grid('on')
plot( timeaxis, xF2a(2,:), LineWidth=1.25 )
plot( timeaxis, xFRa(2,:), 'k--', LineWidth=1.25 )
title('Eccentricity (X) of SC1 and SC2')
legend('SC1','SC2','Reference')

figure(3) % Eccentricity vector Y (unitless)
plot( timeaxis, xF1a(3,:), LineWidth=1.25 )
hold('on')
grid('on')
plot( timeaxis, xF2a(3,:), LineWidth=1.25 )
plot( timeaxis, xFRa(3,:), 'k--', LineWidth=1.25 )
title('Eccentricity (Y) of SC1 and SC2')
legend('SC1','SC2','Reference')

figure(4) % Inclination angle (degrees)
plot( timeaxis, rad2deg( xF1a(4,:)), LineWidth=1.25 )
hold('on')
grid('on')
plot( timeaxis, rad2deg( xF2a(4,:)), LineWidth=1.25 )
plot( timeaxis, rad2deg( xFRa(4,:)), 'k--', LineWidth=1.25 )
title('Inclination (rad) of SC1 and SC2')
legend('SC1','SC2','Reference')

figure(5) % Right angle of ascending node (degrees)
plot( timeaxis, rad2deg( xF1a(5,:)), LineWidth=1.25 )
hold('on')
grid('on')
plot( timeaxis, rad2deg( xF2a(5,:)), LineWidth=1.25 )
plot( timeaxis, rad2deg( xFRa(5,:)), 'k--', LineWidth=1.25 )
title('RAAN (rad) of SC1 and SC2')
legend('SC1','SC2','Reference')

figure(6) % Argument of latitude (degrees)
plot( timeaxis, rad2deg( xF1a(6,:)), LineWidth=1.25 )
hold('on')
grid('on')
plot( timeaxis, rad2deg( xF2a(6,:)), LineWidth=1.25 )
plot( timeaxis, rad2deg( xFRa(6,:)), 'k--', LineWidth=1.25 )
title('Arg of Latitude (rad) of SC1 and SC2')
legend('SC1','SC2','Reference')

figure(7) % Delta-V Usage of SC1/2 at each time step
plot( timeaxis, uF1a', LineWidth=1.25 )
hold('on')
grid('on')
plot( timeaxis, uF2a', LineWidth=1.25 )
title('Delta-V Usage each time step')
legend('SC1-V_R','SC1-V_T','SC1-V_N',...
       'SC2-V_R','SC2-V_T','SC2-V_N')

figure(8) % Cumulative Delta-V Usage of SC1/2
plot( timeaxis, dv1a(2:end), LineWidth=1.25 )
hold('on')
grid('on')
plot( timeaxis, dv2a(2:end), LineWidth=1.25 )
title('Cumulative Delta-V Usage in total')
legend('SC1','SC2')



%% Primary Decentralized MPC Optimization Process (Analytical).
% =========================================================================

function [xf1, xf2, xfR, uf1, uf2] = run_DecMPC( xk1, xk2, xkR, dt, ...
                                                 N, Q, R, P, u_lb, u_ub)

    % Input: xk1, xk2, are 6xN matrices detailing the states of SC1/2
    % Initialize MPC global and local variables.
    u10 = zeros(3,N); % Initialize RTN control ΔV for SC1 [m/s]
    u20 = zeros(3,N); % Initialize RTN control ΔV for SC2 [m/s]
    uf1 = u10;
    uf2 = u20;
    xf1 = xk1(:,1); 
    xf2 = xk2(:,1);
    xfR = xkR(:,1);
    mu = 3.986004415e14;
    Rinv = inv(R);

    % Re-initialize (re-update) the initial state to length 1.
    xf1 = xf1(:,1); % Take just the first vector.
    xf2 = xf2(:,1); % Take just the first vector.
    xfR = xfR(:,1); % Take just the first vector.

    % Populate the states of both spacecraft first.
    for i = 1:N
        xf1(:,i+1) = nonlinear_dynamics( xf1(:,i), uf1(:,i), dt );
        xf2(:,i+1) = nonlinear_dynamics( xf2(:,i), uf2(:,i), dt );
        xfR(:,i+1) = nonlinear_dynamics( xfR(:,i), zeros(3,1), dt );
    end

    % (Re)-populate states based on control from previous iteration.
    for i = 1:N

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
        x1_0_term = -1 * B1' * Q * ( xf1(:,i+1) - xfR(:,i+1) );
        x1_F_term = -1 * B1' * A1stack * P * (xf1(:,end)-xfR(:,end));
        uf1(:,i) = x1_0_term + x1_F_term;
        uf1(:,i) = Rinv * uf1(:,i);
        
        % Compute the optimal control vector for SC2
        x2_0_term = -1 * B2' * Q * ( xf2(:,i+1) - xfR(:,i+1) );
        x2_F_term = -1 * B2' * A2stack * P * (xf2(:,end)-xfR(:,end));
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


