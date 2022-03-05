% Basic simulation parameters
m = 3; % Prediction horizon
dt = 1.0; % Time step

% Initialize Keplerian elements (a,e,i,w,R,nu).
kep1 = [7000000, 0.2, 45, 60, 90, 180];
kep2 = [6995000, 0.22, 47, 62, 91, 175];

% Convert orbit elements for Tomi's dynamics function
%  --> 1. a  - semi-major axis [m]
%  --> 2. ex - x-component of eccentricity vector [-]
%  --> 3. ey - y-component of eccentricity vector [-]
%  --> 4. i  - inclination [rad]
%  --> 5. Om - right ascension of the ascending node [rad]
%  --> 6. u  - argument of latitude [rad]

oe1 = [ kep1(1) ...
        kep1(2)*cosd( kep1(4) ) ...
        kep1(2)*sind( kep1(4) ) ...
        deg2rad( kep1(3) ) ...
        deg2rad( kep1(5) ) ...
        deg2rad( kep1(6) + kep1(4) ) ];

oe2 = [ kep2(1) ...
        kep2(2)*cosd( kep2(4) ) ...
        kep2(2)*sind( kep2(4) ) ...
        deg2rad( kep2(3) ) ...
        deg2rad( kep2(5) ) ...
        deg2rad( kep2(6) + kep2(4) ) ];

% Initialize the state and costs into one vector
x = [ oe1 oe2 zeros()]

% Run optimization here, input initial du, output final du
x = patternsearch(fun,x0)

% Apply the very first du vector, propagate the dynamics
oe1 = nonlinear_dynamics( oe1, u1, dt );

%% Section of the code where some cost function is evaluated
% (Michael is working on this)
function J = dummy_cost(x)
    % Structure of the inputs
    u1 = x(1:)



    % m ------> horizon length
    % u1 -----> mx3 matrix of control inputs for SC1
    % u2 -----> mx3 matrix of control inputs for SC2
    % oe1_i --> vector of [a,ex,ey,i,Om,u] for SC1
    % oe2_i --> vector of [a,ex,ey,i,Om,u] for SC2

    % First, evaluate the cost for SC1
    J = 0.0;
    oe1_l = oe1_i;
    oe2_l = oe2_i;
    for i = 1:m
        J = J + norm( oe1_l - oe2_l );
        oe1_l = nonlinear_dynamics( oe1_l, u1(i,:), dt );
        oe2_l = nonlinear_dynamics( oe2_l, u2(i,:), dt );
    end
    J = J + sum(vecnorm(u1'));
    J = J + sum(vecnorm(u2'));
end