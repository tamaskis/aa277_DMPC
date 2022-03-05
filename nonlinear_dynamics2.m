%==========================================================================
%
% nonlinear_dynamics  Discrete-time nonlinear dynamics of a spacecraft in
% a near-circular orbit subject to impulsive maneuvers.
%
%   x_next = nonlinear_dynamics(x,u)
%
% Author: Tamas Kis
% Last Update: 2022-03-04
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   x       - ((6 or 7)×1 double) state vector at current sample time
%               --> 1. a  - semi-major axis [m]
%               --> 2. ex - x-component of eccentricity vector [-]
%               --> 3. ey - y-component of eccentricity vector [-]
%               --> 4. i  - inclination [rad]
%               --> 5. Om - right ascension of the ascending node [rad]
%               --> 6. u  - argument of latitude [rad]
%               --> 7. dV - (OPTIONAL) argument of latitude [rad]
%   u       - (3×1 double) control input at current sample time
%               --> 1. dvR - impulsive maneuver in radial direction [m/s]
%               --> 2. dvT - impulsive maneuver in tang. direction [m/s]
%               --> 3. dvN - impulsive maneuver in normal direction [m/s]
%   dt      - (1×1 double) time step [s]
%
% -------
% OUTPUT:
% -------
%   x_next  - (7×1 double) state vector at next sample time
%
%==========================================================================
function x_next = nonlinear_dynamics(x,u,dt)
    
    % unpacks control input
    dvR = u(1);
    dvT = u(2);
    dvN = u(3);

    % unpacks state vector (note: name "u" assigned to new variable)
    a = x(1);
    ex = x(2);
    ey = x(3);
    i = x(4);
    Om = x(5);
    u = x(6);
    if (length(x) == 7), dV = x(7); end

    % Earth gravitational parameter [m^3/s^2]
    mu = 3.986004415e14;

    % mean motion [rad/s]
    n = sqrt(mu/(a)^3);
    
    % changes in orbital elements due to impulse maneuver
    da = (2/n)*dvT;
    nF = sqrt(mu/(a+da)^3); % Updated n
    dex = (sin(u)/(n*a))*dvR+(2*cos(u)/(n*a))*dvT;
    dey = (-cos(u)/(n*a))*dvR+(2*sin(u)/(n*a))*dvT;
    di = (cos(u)/(n*a))*dvN;
    dOm = (sin(u)/(n*a*sin(i)))*dvN;
    du = nF*dt-(2/(n*a))*dvR-((sin(u)*cot(i))/(n*a))*dvN;

    % orbital elements at next sample time
    a_next = a+da;
    ex_next = ex+dex;
    ey_next = ey+dey;
    i_next = i+di;
    Om_next = Om+dOm;
    u_next = u+du;

    % cumulative expended dV at next sample time
    dV_next = dV+dvR+dvT+dvN;

    % packages state vector at next sample time
    x_next = [a_next;
              ex_next;
              ey_next;
              i_next;
              Om_next;
              u_next];

    % adds dV_next to state vector if state vector was input with dV
    if (length(x) == 7), x_next = [x_next;dV_next]; end

end