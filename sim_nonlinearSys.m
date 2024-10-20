% clear variables
% close all
% clc

SIMULATION_TIME = 20;
STEP_SIZE = 5e-3; % step size
%% Systen Parameters
identified_para = matfile("data/identified_parameters.mat");
a1 = identified_para.a1;
a2 = identified_para.a2;
b1 = identified_para.b1;
b2 = identified_para.b2;

g  = 9.81;
L1 = 9.20e-02;  % 9.2cm

%% Initial Condition
% x0 = [-0.2 0 0 0]'; % comment out this line to run this script with different x0
% x0 = [-0.2693 1.4963 0.2051 -0.8728]';
% K = [10 1 10 1];
K = [5.9934 0.5486 6.8171 0.3307]; % K_1 with R=100 nonlinear iteration ppt:20241016
K = [4.5409 0.2752 5.4382 0.0631]; % K_2 with R=100 nonlinear iteration ppt:20241016
K = [3.9809 0.1691 4.9928 -0.0086]; % K_3 with R=100 nonlinear iteration ppt:20241016
K = [3.8471 0.1344 4.9866 0.0033]; % K_4 with R=100 nonlinear iteration ppt:20241016

sim_NL = sim('sim_nonlinear_pend.slx');

%% Plot
legend_name_NL = ["motor angle sim NL","motor speed sim NL","pendulum angle sim NL","pendulum speed sim NL"];

figure("Name","Nonlinear Simulation")
sgtitle("Nonlinear Simulation",'Interpreter','latex')
generatePhyiscalDataPlot("a",sim_NL.state4dim.Data',sim_NL.tsim',legend_name_NL)
generateStateDataPlot("Simulation State Data",sim_NL.state4dim.Data',sim_NL.tsim')

% ---- compare with experiment data
% hold on
% generatePhyiscalDataPlot("a",sim_NL.state4dim.Data',sim_NL.tsim,legend_name_NL)
% ----

% %% Select Data from Linear Sim Data
% X = cell(0,0);            % comment out this line to add data
% T = cell(0,0);            % comment out this line to add data
% % load("data\data_temp.mat")  % comment out this line when add the first data group
% 
% % select time interval and segment number and segment offset size
% dataPros_SimNL2Datadriven = dataProcessing(sim_NL.state4dim.Data',sim_NL.tsim',[]);
% T_int = [0,0.05];
% offset_size = 0.05;
% segment_number = 10;
% for i = 1:segment_number
%     [x_dd, t_dd, ~] = dataPros_SimNL2Datadriven.getTimeIntervalData(T_int);
%     X = [X; x_dd];
%     T = [T; t_dd];
%     T_int = T_int + offset_size;
% end
% 
% %% Plot Data-driven Data
% % figure("Name","All Data-driven Data")
% % sgtitle("All Data-driven Data - NL System",'Interpreter','latex')
% for i = 1:size(X,1)
%     subplot(size(X,1)/segment_number,segment_number,i)
%     hold on
%     plot(T{i},X{i},'rsquare-',MarkerSize=6)
%     xlabel('$t$','Interpreter','latex')
%     ylabel('$x$','Interpreter','latex')
%     xlim([T{i}(1),T{i}(end)])
%     grid on
% end
% 
% %% Data-driven Solution
% E = [1 0 0 0;
%      0 1 0 0;
%      0 0 1 0;
%      0 -L1 0 a2 ];
% F = [ 0  1  0  0
%       0 -a1 0  0
%       0  0  0  1
%       0  0 -g -b2 ];
% G = [ 0
%       b1
%       0
%       0 ];
% A = E \ F;
% B = E \ G;
% q1 = 100;
% q3 = 1000;
% Q = diag([q1 0 q3 0]);
% R = 1;
% sys = linearSys(A,B,x0,Q,R);
% 
% l = 0;
% dd_sim_NL = ddLyap(X, l, STEP_SIZE, K, sys); % solve data driven Lyapunov equation
% Pi = dd_sim_NL.Pi;
% Kip1 = dd_sim_NL.Kip1;
% data_end_point = dd_sim_NL.data_matrix_end_point;
% data_int = dd_sim_NL.data_array_integral;
% 
% % calculate error matrix
% Delta_i = Pi - dd_sim_NL.Pi_lyap;
% norm_Delta_i = norm(Delta_i);
% fprintf(['|Delta_i| = ', '\n'])
% fprintf([num2str(norm(Delta_i)), '\n'])


%% --------------------------- Functions --------------------------
function generatePhyiscalDataPlot(title,x,t,LegendName)
% figure("Name",title)
% sgtitle(title,'Interpreter','latex')
subplot(3,1,1)
hold on
plot(t,x(1,:),'--','DisplayName',LegendName(1),LineWidth=2)
plot(t,x(3,:),'--','DisplayName',LegendName(3),LineWidth=2)
legend('Interpreter','latex',Location="best")
xlim([t(1),t(end)])
grid on

subplot(3,1,2)
hold on
plot(t,x(2,:),'--','DisplayName',LegendName(2),LineWidth=2)
plot(t,x(4,:),'--','DisplayName',LegendName(4),LineWidth=2)
legend('Interpreter','latex',Location="best")
xlim([t(1),t(end)])
grid on
end

function generateStateDataPlot(title,x,t)
figure("Name",title)
sgtitle(title,'Interpreter','latex')
plot(t,x,'.-',MarkerSize=10);
legend('$\theta_1$','$\dot \theta_1$','$\phi_2$','$\dot \phi_2$','Interpreter','latex')
xlabel('$t$','Interpreter','latex')
ylabel('$x$','Interpreter','latex')
xlim([t(1),t(end)])
grid on
ax = gca;
ax.FontSize = 14;
end

