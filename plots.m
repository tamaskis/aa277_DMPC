%% Plots for DMPC Project
% =========================================================================
% AA277  |  Luke Neise, Samuel Low, Michael Ying, Tamas Kis

clc; clear all; close all;

% load data
DMPC_analytical = struct2array(load('data/DMPC_analytical'));
DecMPC_analytical = struct2array(load('data/DecMPC_analytical'));
DMPC_numerical = struct2array(load('data/DMPC_numerical'));
DecMPC_numerical = struct2array(load('data/DecMPC_numerical'));

% load plot parameters
pp = PLOT_PARAMETERS;



%% PLOTS

% ------------------------------
% Plots for analytical solution.
% ------------------------------

%plot_result(DMPC_analytical,DecMPC_analytical,pp,'Analytical',24)

% -----------------------------
% Plots for numerical solution.
% -----------------------------

plot_result(DMPC_numerical,DecMPC_numerical,pp,'Numerical',24)

% -------------------------------------------------------------
% Representative plot of delta-V usage at each time step [m/s].
% -------------------------------------------------------------

% figure('Position',[540,100,700,300]);
% hold on;
% plot(DMPC_numerical.t,DMPC_numerical.ctrl1(1,:),'LineWidth',1.5,...
%     'Color',pp.matlab_blue);
% plot(DMPC_numerical.t,DMPC_numerical.ctrl2(1,:),'LineWidth',1.5,...
%     'Color',pp.matlab_blue,'LineStyle','--');
% plot(DMPC_numerical.t,DMPC_numerical.ctrl1(2,:),'LineWidth',1.5,...
%     'Color',pp.matlab_red);
% plot(DMPC_numerical.t,DMPC_numerical.ctrl2(2,:),'LineWidth',1.5,...
%     'Color',pp.matlab_red,'LineStyle','--');
% plot(DMPC_numerical.t,DMPC_numerical.ctrl1(3,:),'LineWidth',1.5,...
%     'Color',pp.matlab_yellow);
% plot(DMPC_numerical.t,DMPC_numerical.ctrl2(3,:),'LineWidth',1.5,...
%     'Color',pp.matlab_yellow,'LineStyle','--');
% hold off;
% grid on;
% xlim([0,0.75]);
% ylabel('$\Delta V\;[\mathrm{m/s}]$ (at each time step)',...
%     'Interpreter','latex','FontSize',18);
% title('\boldmath$\Delta V$ \textbf{Usage}','Interpreter','latex',...
%     'FontSize',18);
% legend('$(\Delta v_{R})_{1}$','$(\Delta v_{R})_{2}$',...
%     '$(\Delta v_{T})_{1}$','$(\Delta v_{T})_{2}$',...
%     '$(\Delta v_{N})_{1}$','$(\Delta v_{N})_{2}$','Interpreter',...
%     'latex','FontSize',14,'Location','southeast');



%% FUNCTION TO PLOT RESULTS

function plot_result(DMPC,DecMPC,pp,an_or_num,tmax)

    % initialize figure
    figure('Position',[540,100,700,800]);
    
    % semi-major axis [km]
    subplot(6,1,1);
    hold on;
    plot(DMPC.t,DMPC.a1,'LineWidth',1.5,'Color',pp.matlab_blue)
    plot(DecMPC.t,DecMPC.a1,'LineWidth',1.5,'Color',pp.matlab_light_blue);
    plot(DMPC.t,DMPC.a2,'LineWidth',1.5,'Color',pp.matlab_red)
    plot(DecMPC.t,DecMPC.a2,'LineWidth',1.5,'Color',pp.matlab_light_red);
    hold off;
    grid on;
    xlim([0,tmax]);
    ylabel('$a\;[\mathrm{km}]$','Interpreter','latex','FontSize',18);
    title({"\textbf{Orbital Element Trajectories ("+an_or_num+...
        " Solution)}",''},'Interpreter','latex','FontSize',18);
    legend('spacecraft 1 (DMPC)','spacecraft 1 (DecMPC)',...
        'spacecraft 2 (DMPC)','spacecraft 2 (DecMPC)','Interpreter',...
        'latex','FontSize',12,'Position',[0.28,0.41,1,1]);

    % x-component of eccentricity vector [-]
    subplot(6,1,2);
    hold on;
    plot(DMPC.t,DMPC.ex1,'LineWidth',1.5,'Color',pp.matlab_blue)
    plot(DecMPC.t,DecMPC.ex1,'LineWidth',1.5,'Color',pp.matlab_light_blue);
    plot(DMPC.t,DMPC.ex2,'LineWidth',1.5,'Color',pp.matlab_red)
    plot(DecMPC.t,DecMPC.ex2,'LineWidth',1.5,'Color',pp.matlab_light_red);
    hold off;
    grid on;
    xlim([0,tmax]);
    ylabel('$e_{x}$','Interpreter','latex','FontSize',18);

    % y-component of eccentricity vector [-]
    subplot(6,1,3);
    hold on;
    plot(DMPC.t,DMPC.ey1,'LineWidth',1.5,'Color',pp.matlab_blue)
    plot(DecMPC.t,DecMPC.ey1,'LineWidth',1.5,'Color',pp.matlab_light_blue);
    plot(DMPC.t,DMPC.ey2,'LineWidth',1.5,'Color',pp.matlab_red)
    plot(DecMPC.t,DecMPC.ey2,'LineWidth',1.5,'Color',pp.matlab_light_red);
    hold off;
    grid on;
    xlim([0,tmax]);
    ylabel('$e_{y}$','Interpreter','latex','FontSize',18);
    
    % inclination [deg]
    subplot(6,1,4);
    hold on;
    plot(DMPC.t,DMPC.i1,'LineWidth',1.5,'Color',pp.matlab_blue);
    plot(DecMPC.t,DecMPC.i1,'LineWidth',1.5,'Color',pp.matlab_light_blue);
    plot(DMPC.t,DMPC.i2,'LineWidth',1.5,'Color',pp.matlab_red);
    plot(DecMPC.t,DecMPC.i2,'LineWidth',1.5,'Color',pp.matlab_light_red);
    hold off;
    grid on;
    xlim([0,tmax]);
    ylabel('$i\;[{}^{\circ}]$','Interpreter','latex','FontSize',18);
    
    % RAAN [deg]
    subplot(6,1,5);
    hold on;
    plot(DMPC.t,DMPC.Om1,'LineWidth',1.5,'Color',pp.matlab_blue);
    plot(DecMPC.t,DecMPC.Om1,'LineWidth',1.5,'Color',pp.matlab_light_blue);
    plot(DMPC.t,DMPC.Om2,'LineWidth',1.5,'Color',pp.matlab_red);
    plot(DecMPC.t,DecMPC.Om2,'LineWidth',1.5,'Color',pp.matlab_light_red);
    hold off;
    grid on;
    xlim([0,tmax]);
    ylabel('$\Omega\;[{}^{\circ}]$','Interpreter','latex','FontSize',18);

    % argument of latitude [deg]
    subplot(6,1,6);
    hold on;
    plot(DMPC.t,DMPC.u2-DMPC.u1,'LineWidth',1.5,'Color',pp.matlab_blue);
    plot(DecMPC.t,DecMPC.u2-DecMPC.u1,'LineWidth',1.5,'Color',...
        pp.matlab_light_blue);
    hold off;
    grid on;
    xlim([0,tmax]);
    xlabel('Time $[\mathrm{h}]$','Interpreter','latex','FontSize',18);
    ylabel('$u_{2}-u_{1}\;[{}^{\circ}]$','Interpreter','latex',...
        'FontSize',18);
    legend('DMPC','DecMPC','Interpreter','latex','FontSize',12,...
        'Location','southeast');
    
    % -------------------------------
    % Cumulative delta-V usage [m/s].
    % -------------------------------
    
    figure('Position',[540,100,700,300]);
    hold on;
    plot(DMPC.t,DMPC.dV1,'LineWidth',1.5,'Color',pp.matlab_blue);
    plot(DecMPC.t,DecMPC.dV1,'LineWidth',1.5,'Color',pp.matlab_light_blue);
    plot(DMPC.t,DMPC.dV2,'LineWidth',1.5,'Color',pp.matlab_red);
    plot(DecMPC.t,DecMPC.dV2,'LineWidth',1.5,'Color',pp.matlab_light_red);
    hold off;
    grid on;
    xlim([0,tmax]);
    ylabel('$\Delta V$ (cumulative) $[\mathrm{m/s}]$','Interpreter',...
        'latex','FontSize',18);
    title("\textbf{Cumulative} \boldmath$\Delta V$ \textbf{Usage ("+...
        an_or_num+" Solution)}",'Interpreter','latex','FontSize',18);
    legend('spacecraft 1 (DMPC)','spacecraft 1 (DecMPC)',...
        'spacecraft 2 (DMPC)','spacecraft 2 (DecMPC)','Interpreter',...
        'latex','FontSize',12,'Location','best');
    
    % RTN Position Plot
    figure;
    hold on;
    plot(DMPC.t,DMPC.RTN_2wrt1(1,:),'LineWidth',1.5,'Color',...
        pp.matlab_blue);
    plot(DecMPC.t,DecMPC.RTN_2wrt1(1,:),'LineWidth',1.5,'Color',...
        pp.matlab_light_blue);
    plot(DMPC.t,DMPC.RTN_2wrt1(2,:),'LineWidth',1.5,'Color',...
        pp.matlab_red);
    plot(DecMPC.t,DecMPC.RTN_2wrt1(2,:),'LineWidth',1.5,'Color',...
        pp.matlab_light_red);
    plot(DMPC.t,DMPC.RTN_2wrt1(3,:),'LineWidth',1.5,'Color',...
        pp.matlab_yellow);
    plot(DecMPC.t,DecMPC.RTN_2wrt1(3,:),'LineWidth',1.5,'Color',...
        pp.matlab_light_yellow);
    hold off;
    grid on;
    xlim([0,tmax]);
    xlabel('Time $[\mathrm{h}]$','Interpreter','latex','FontSize',18);
    ylabel('Relative Position $[\mathrm{km}]$','Interpreter','latex',...
        'FontSize',18);
    title("\textbf{Relative RTN Position of Spacecraft 2 w.r.t. Space"+...
        "craft 1}",'Interpreter','latex','FontSize',18);
    legend('$R$ (DMPC)','$R$ (DecMPC)','$T$ (DMPC)','$T$ (DecMPC)',...
        '$N$ (DMPC)','$N$ (DecMPC)','Interpreter','latex','FontSize',12,...
        'Location','best');

end