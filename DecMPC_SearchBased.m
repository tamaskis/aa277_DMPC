%% Decentralized MPC Algorithm for Spacecraft Rendezvous
% AA277  |  Luke Neise, Samuel Low, Michael Ying, Tamas Kis

clear; clc; close all;

%% Exterior Online Simulation

% Set up our spacecraft and reference
x0_sc1 = [6700*1000, 0.002, 0.002, 45.05*pi/180, 45.05*pi/180, 1*pi/180]'; % define initial oe []
x0_sc2 = [6750*1000, -0.002, -0.002, 45*pi/180, 45*pi/180, -1*pi/180]'; % define initial oe []
x0_ref = [6.7261e6, 0.0020273, 0.0010998, 0.786, 0.78634, 0]'; % define reference oe []
n = sqrt(3.986004415e14/x0_ref(1)^3); % mean motion of reference

% Set up MPC timing/horizon
dt = 600; thorizon = [0:dt:1.5*3600]; N = 3; 
dt = 60; thorizon = [0:dt:86400]; N = 5; 

% Simulate the rendezvous
x_sc1 = []; u_sc1 = []; u_sc1_tot = [];
x_sc2 = []; u_sc2 = []; u_sc2_tot = [];
t_runtime = zeros(1,length(thorizon));
t_runtimesum = zeros(1,length(thorizon));
for jj = 1:length(thorizon)
    % Find time-dependent reference
    Xref = repmat(x0_ref,[1,N]);
    for ii = 1:N
        Xref(6,ii) = x0_ref(6) + n*dt*(jj-1+ii-1);
    end

    % Find optimal MPC control
    tic;
    U_opt_sc1 = MPCprocess(x0_sc1,Xref,dt,N);
    U_opt_sc2 = MPCprocess(x0_sc2,Xref,dt,N);
    u1 = U_opt_sc1(1:3,1);
    u2 = U_opt_sc2(1:3,1);
    t_runtime(jj) = toc;
    if jj > 1
        t_runtimesum(jj) = t_runtimesum(jj-1) + t_runtime(jj);
    else
        t_runtimesum(jj) = t_runtime(jj);
    end
    
    % Take step with control
    x1_sc1 = nonlinear_dynamics(x0_sc1,u1,dt);
    x1_sc2 = nonlinear_dynamics(x0_sc2,u2,dt);
    
    % Record
    x_sc1 = [x_sc1, x0_sc1];
    x_sc2 = [x_sc2, x0_sc2];
    u_sc1 = [u_sc1, u1];
    u_sc2 = [u_sc2, u2];
    if jj>1
        u_sc1_tot(jj) = u_sc1_tot(jj-1) + norm(u1);
        u_sc2_tot(jj) = u_sc2_tot(jj-1) + norm(u2);
    else
        u_sc1_tot(jj) = norm(u1);
        u_sc2_tot(jj) = norm(u2);
    end
    
    % Display
    %disp(u1);
    %disp(u2);
    %disp(x0_sc1-x0_sc2);

    % Update
    x0_sc1 = x1_sc1;
    x0_sc2 = x1_sc2;   
    
end

duration = 86400;

timeaxis = linspace( 0, duration, (duration/dt)+1 );

% converts time to hours
timeaxis = timeaxis/3600;

% converts semi-major axes to kilometers
x_sc1(1,:) = x_sc1(1,:)/1000;
x_sc2(1,:) = x_sc2(1,:)/1000;


%%

% Get RTN
RTN_sc2_wrtsc1 = oe2rtnback(x_sc1,x_sc2);
RTN_sc2_wrtsc1 = RTN_sc2_wrtsc1.';

% converts RTN to km
RTN_sc2_wrtsc1 = RTN_sc2_wrtsc1/1000;


%% Saving data

% time [h]
t = timeaxis;

% semi-major axes [km]
a1 = x_sc1(1,:);
a2 = x_sc2(1,:);

% x-component of eccentricity vector [-]
ex1 = x_sc1(2,:);
ex2 = x_sc2(2,:);

% y-component of eccentricity vector [-]
ey1 = x_sc1(3,:);
ey2 = x_sc2(3,:);

% inclination [deg]
i1 = x_sc1(4,:);
i2 = x_sc2(4,:);

% RAAN [deg]
Om1 = x_sc1(5,:);
Om2 = x_sc2(5,:);

% argument of latitude [deg]
u1 = x_sc1(6,:);
u2 = x_sc2(6,:);

% control effort [m/s]
ctrl1 = u_sc1;
ctrl2 = u_sc2;

% cumulative Delta-V [m/s]
dV1 = u_sc1_tot;
dV2 = u_sc2_tot;

% RTN positions [km]
RTN_2wrt1 = RTN_sc2_wrtsc1;

% trims all arrays so they are same length as for analytical case (only
% differs by 1 element) and packages them into a structure
DecMPC_numerical.t = t(1:1440);
DecMPC_numerical.a1 = a1(1:1440);
DecMPC_numerical.a2 = a2(1:1440);
DecMPC_numerical.ex1 = ex1(1:1440);
DecMPC_numerical.ex2 = ex2(1:1440);
DecMPC_numerical.ey1 = ey1(1:1440);
DecMPC_numerical.ey2 = ey2(1:1440);
DecMPC_numerical.i1 = i1(1:1440);
DecMPC_numerical.i2 = i2(1:1440);
DecMPC_numerical.Om1 = Om1(1:1440);
DecMPC_numerical.Om2 = Om2(1:1440);
DecMPC_numerical.u1 = u1(1:1440);
DecMPC_numerical.u2 = u2(1:1440);
DecMPC_numerical.ctrl1 = ctrl1(:,1:1440);
DecMPC_numerical.ctrl2 = ctrl2(:,1:1440);
DecMPC_numerical.dV1 = dV1(1:1440);
DecMPC_numerical.dV2 = dV2(1:1440);
DecMPC_numerical.RTN_2wrt1 = RTN_2wrt1(:,1:1440);

% saves simulation data
save('data/DecMPC_numerical.mat','DecMPC_numerical');


%% Plots
% RTN convergence from reference of final orbit
fig1 = figure(1);
plot(thorizon,RTN_sc2_wrtsc1(:,1));
hold on;
plot(thorizon,RTN_sc2_wrtsc1(:,2));
plot(thorizon,RTN_sc2_wrtsc1(:,3));
xlabel('Simulation Time [s]');
ylabel('Distance [m]');
title('Spacecraft 2 Distance from Spacecraft 1 in RTN');
legend('R','T','N');

% Control effort over time
fig2 = figure(2);
hold on;
plot(thorizon,u_sc1(1,:),'k','Linestyle','none','Marker','o');
plot(thorizon,u_sc1(2,:),'k','Linestyle','none','Marker','*');
plot(thorizon,u_sc1(3,:),'k','Linestyle','none','Marker','d');
plot(thorizon,u_sc2(1,:),'b','Linestyle','none','Marker','o');
plot(thorizon,u_sc2(2,:),'b','Linestyle','none','Marker','*');
plot(thorizon,u_sc2(3,:),'b','Linestyle','none','Marker','d');
title('Decentralized MPC Control Effort at Each Step');
xlabel('Simulation Time [s]');
ylabel('Control Effort [m/s]');
legend('sc1 dvR','sc1 dvT','sc1 dvN','sc2 dvR','sc2 dvT','sc2 dvN','location','East');

fig3 = figure(3);
hold on;
plot(thorizon,u_sc1_tot,'k--');
plot(thorizon,u_sc2_tot,'b--');
title('Cumulative Decentralized MPC Control Effort');
xlabel('Simulation Time [s]');
ylabel('Control Effort [m/s]');
legend('sc1 cumulative dv','sc2 cumulative dv','location','East');

% Runtime
fig4 = figure(4);
hold on;
plot(thorizon,t_runtime,'k','Linestyle','none','Marker','o');
title('Decentralized MPC Runtime at Each Step');
xlabel('Simulation Time [s]');
ylabel('Step Runtime [s]');
legendstr = sprintf('Total Runtime - %.1f sec',t_runtimesum(end));
legend(legendstr);

% oe evolution
figure;
hold on;
plot(thorizon,x_sc1(1,:));
plot(thorizon,x_sc2(1,:));
title('Semimajor Axis vs. Time');
xlabel('Simulation Time [s]');
ylabel('Semimajor Axis [m]');
legend('sc1','sc2');

figure;
hold on;
plot(thorizon,x_sc1(2,:));
plot(thorizon,x_sc2(2,:));
title('ey vs. Time');
xlabel('Simulation Time [s]');
ylabel('ey');
legend('sc1','sc2');

figure;
hold on;
plot(thorizon,x_sc1(3,:));
plot(thorizon,x_sc2(3,:));
title('ex vs. Time');
xlabel('Simulation Time [s]');
ylabel('ex');
legend('sc1','sc2');

figure;
hold on;
plot(thorizon,x_sc1(4,:));
plot(thorizon,x_sc2(4,:));
title('Inclination vs. Time');
xlabel('Simulation Time [s]');
ylabel('Inclination [rad]');
legend('sc1','sc2');

figure;
hold on;
plot(thorizon,x_sc1(5,:));
plot(thorizon,x_sc2(5,:));
title('RAAN vs. Time');
xlabel('Simulation Time [s]');
ylabel('RAAN [rad]');
legend('sc1','sc2');

figure;
hold on;
plot(thorizon,x_sc1(6,:));
plot(thorizon,x_sc2(6,:));
title('Argument of Latitude vs. Time');
xlabel('Simulation Time [s]');
ylabel('Argument of Latitude [rad]');
legend('sc1','sc2');


%% MPC Process and Supplement Functions

% Run Model Predictive Control top level
function U_opt = MPCprocess(x0,Xref,dt,N)

% Initialize MPC details
% Q = eye(6); % define state error cost matrix
% R = 0.1*eye(3); % define control effort cost matrix
% Qf = eye(6); % define terminal state error cost matrix 
Q = 0.1*eye(6);
R = 0.001*eye(3);
Qf = (N^2)*eye(6);
Q(1,1) = 0.2*N / (0.0175 * (x0(1)+Xref(1,1))/2);
Qf(1,1) = 0.2*N / (0.0175 * (x0(1)+Xref(1,1))/2);
u_lb = -3*ones(3*N,1); % define deltaV lower bound [m/s]
u_ub = 3*ones(3*N,1); % define deltaV upper bound [m/s]
U0 = ones(3*N,1); % define RTN control deltaVs [m/s]

% Run MPC
options = optimset('Display','off');
fun = @(U) MPCsim(x0,Xref,U,Q,R,Qf,N,dt);
U_opt = fmincon(fun,U0,[],[],[],[],u_lb,u_ub,[],options);

end

% Simulate horizon cost
function totalcost = MPCsim(x0,Xref,U,Q,R,Qf,N,dt)

% Prep storage
xi_vec = [x0];
totalcost = 0;

% Predict horizon
for ii = 1:N-1
    ui = U((3*ii - 2):3*ii);
    xicurr = xi_vec(:,ii);
    xiref = Xref(:,ii);
    xinext = nonlinear_dynamics(xicurr,ui,dt);
    xi_vec = [xi_vec, xinext];
    subcost = cost(xiref,xicurr,ui,Q,R);
    totalcost = totalcost + subcost;
end
xicurr = xi_vec(:,N);
xiref = Xref(:,N);
subcost = termcost(xiref,xicurr,Qf);
totalcost = totalcost + subcost;

end

% Compute next state
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

% Compute cost of state/control
function J = cost(xref, x, u, Q, R)
dx = xref - x;
J = dx'*Q*dx + u'*R*u;
end

% Compute terminal cost
function J = termcost(xref, x, Qf)
dx = xref - x;
J = dx'*Qf*dx;
end

% Compute RTN position
function RTN_sc2_wrtsc1 = oe2rtnback(x_sc1,x_sc2)

% Go from chief and deputy to RTN delta
RTN_sc2_wrtsc1 = zeros(length(x_sc1),3);
for ii = 1:length(x_sc1)
    %dx_sc2_wrtsc1 = elements_to_RTN_Luke(x_sc1(:,ii),x_sc2(:,ii));
    dx_sc2_wrtsc1 = elements_to_RTN(x_sc1(:,ii),x_sc2(:,ii));
    RTN_sc2_wrtsc1(ii,:) = dx_sc2_wrtsc1(1:3);
end

end

% Go orbital elements to RTN
function dx = elements_to_RTN_Luke(oe_sc1,oe_sc2)

% Extract
a1 = oe_sc1(1); ex1 = oe_sc1(2); ey1 = oe_sc1(3); i1 = oe_sc1(4); RAAN1 = oe_sc1(5); u1 = oe_sc1(6);
a2 = oe_sc2(1); ex2 = oe_sc2(2); ey2 = oe_sc2(3); i2 = oe_sc2(4); RAAN2 = oe_sc2(5); u2 = oe_sc2(6);
e1 = sqrt(ex1^2 + ey1^2); e2 = sqrt(ex2^2 + ey2^2); w1 = atan2(ey1,ex1); w2 = atan2(ey2,ex2);
w1 = wrapTo2Pi(w1);
w2 = wrapTo2Pi(w2);
v1 = u1-w1; v2 = u2-w2;
v1 = wrapTo2Pi(v1);
v2 = wrapTo2Pi(v2);
E1 = v1; E2 = v2;

% Go to RTN
[reci1,veci1] = oe2eci(a1,e1,i1,RAAN1,w1,E1);
[reci2,veci2] = oe2eci(a2,e2,i2,RAAN2,w2,E2);
R_eci2rtn1 = eci2rtn_matrix(reci1,veci1);
R_eci2rtn2 = eci2rtn_matrix(reci2,veci2);
rrtn1 = R_eci2rtn1*reci1;
rrtn2 = R_eci2rtn1*reci2;
dx = rrtn2-rrtn1;

end