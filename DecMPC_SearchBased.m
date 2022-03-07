%% Decentralized MPC Algorithm for Spacecraft Rendezvous
% AA277  |  Luke Neise, Samuel Low, Michael Ying, Tamas Kis


%% Exterior Online Simulation

% Set up our spacecraft and reference
x0_sc1 = [6700*1000, 0, 0, 45*pi/180, 0, 2*pi/180]'; % define initial oe []
x0_sc2 = [6750*1000, 0, 0, 45*pi/180, 0, -2*pi/180]'; % define initial oe []
x0_ref = [6.7261e6, 0.00035547, 0.002252, 0.7854, 8.669e-09, 0]'; % define reference oe []
n = sqrt(3.986004415e14/x0_ref(1)^3); % mean motion of reference

% Set up MPC timing/horizon
dt = 600; thorizon = [0:dt:1.5*3600]; N = 3; 

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
    disp(u1);
    disp(u2);
%     disp(x0_sc1-Xref(:,1));
%     disp(x0_sc2-Xref(:,1));
    disp(x0_sc1-x0_sc2);

    % Update
    x0_sc1 = x1_sc1;
    x0_sc2 = x1_sc2;   
    
end

% Plots

% RTN convergence from reference of final orbit


% Control effort over time
fig2 = figure(2);
hold on;
plot(thorizon,u_sc1(1,:),'k','Linestyle','none','Marker','o');
plot(thorizon,u_sc1(2,:),'k','Linestyle','none','Marker','*');
plot(thorizon,u_sc1(3,:),'k','Linestyle','none','Marker','d');
plot(thorizon,u_sc2(1,:),'b','Linestyle','none','Marker','o');
plot(thorizon,u_sc2(2,:),'b','Linestyle','none','Marker','*');
plot(thorizon,u_sc2(3,:),'b','Linestyle','none','Marker','d');
plot(thorizon,u_sc1_tot,'k--');
plot(thorizon,u_sc2_tot,'b--');
title('Decentralized MPC Control Effort');
xlabel('Simulation Time [s]');
ylabel('Control Effort [m/s]');
legend('sc1 dvR','sc1 dvT','sc1 dvN','sc2 dvR','sc2 dvT','sc2 dvN','sc1 cumulative dv','sc2 cumulative dv','location','East');

% Initial and final orbits


% Runtime
fig4 = figure(4);
hold on;
plot(thorizon,t_runtime,'k','Linestyle','none','Marker','o');
plot(thorizon,t_runtimesum,'k--');
title('Decentralized MPC Runtime at Each Step');
xlabel('Simulation Time [s]');
ylabel('Step Runtime [s]');
legend('Step Runtime','Cumulative Runtime','location','East');


%% MPC Process and Supplement Functions

% Run Model Predictive Control top level
function U_opt = MPCprocess(x0,Xref,dt,N)

% Initialize MPC details
Q = eye(6); % define state error cost matrix
R = 0.1*eye(3); % define control effort cost matrix
Qf = eye(6); % define terminal state error cost matrix 
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