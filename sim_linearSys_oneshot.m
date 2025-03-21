% clear variables
% close all
% clc

%% Define system

identified_para = matfile("data/identified_parameters.mat");
a1 = identified_para.a1;
a2 = identified_para.a2;
b1 = identified_para.b1;
b2 = identified_para.b2;

g  = 9.81;
L1 = 9.20e-02;  % 9.2cm

E = [1 0 0 0;
     0 1 0 0;
     0 0 1 0;
     0 -L1 0 a2 ];
F = [ 0  1  0  0
      0 -a1 0  0
      0  0  0  1
      0  0 -g -b2 ];
G = [ 0
      b1
      0
      0 ];
A = E \ F;
B = E \ G;
% B = [1 1 1 1]';
nx = size(A,2); %state dimension
nu = size(B,2); %input dimension
q1 = 100;
q3 = 1000;
Q = diag([q1 0 q3 0]);
R = 100;

% Generate Raw Simulation Data
STEP_SIZE = 5e-3; % step size

% IC
% x0 = [-0.2 0 0.2 0]'; % comment out this line to run this script with different x0

% FB Gain
% K = [10 1 10 1];
% K = [10 1 1 -0.1];
% K = [5.9934 0.5486 6.8171 0.3307]; % K_1 with R=100 nonlinear iteration ppt:20241016
% K = [4.5409 0.2752 5.4382 0.0631]; % K_2 with R=100 nonlinear iteration ppt:20241016
% K = [3.9809 0.1691 4.9928 -0.0086]; % K_3 with R=100 nonlinear iteration ppt:20241016
% K = [3.8471 0.1344 4.9866 0.0033]; % K_4 with R=100 nonlinear iteration ppt:20241016

if ~all(eig(A - B * K)<0)
    fprintf("K0 is NOT a stabilizing gain! \n")
    return
end
sys = linearSys(A,B,x0,Q,R);

% Generate sim data
T_sim = 10;
[x_sim_linear,u_sim_linear,t_sim_linear] = sys.ResponseFromGain(K,STEP_SIZE,T_sim);
% generateStateDataPlot("Simulation State Data",x_sim_raw,t_sim_raw)

% ------ use this to compare with experiment data
% legend_name_linear = ["motor angle sim linear","motor speed sim linear","pendulum angle sim linear","pendulum speed sim linear"];
% hold on
% generatePhyiscalDataPlot("a",x_sim_linear,t_sim_linear,legend_name_linear) 
% ------

%% Select Data from Linear Sim Data
% X = cell(0,0);            % comment out this line to add data
% T = cell(0,0);            % comment out this line to add data
% % load("data\data_temp.mat")  % comment out this line when add the first data group
% 
% % select time interval and segment number and segment offset size
% dataPros_SimLinear2Datadriven = dataProcessing(x_sim_linear,t_sim_linear,[]);
% T_int = [0,0.05];
% offset_size = 0.05;
% segment_number = 10;
% for i = 1:segment_number
%     [x_dd, t_dd, ~] = dataPros_SimLinear2Datadriven.getTimeIntervalData(T_int);
%     X = [X; x_dd];
%     T = [T; t_dd];
%     T_int = T_int + offset_size;
% end
% 
% %% Plot Data-driven Data
% figure("Name","All Data-driven Data")
% sgtitle("All Data-driven Data - Linear System",'Interpreter','latex')
% for i = 1:size(X,1)
%     subplot(size(X,1)/segment_number,segment_number,i)
%     plot(T{i},X{i},'b.-',MarkerSize=10)
%     xlabel('$t$','Interpreter','latex')
%     ylabel('$x$','Interpreter','latex')
%     xlim([T{i}(1),T{i}(end)])
%     grid on
% end
% %% Data-driven Solution
% l = 0;
% dd_sim_linear = ddLyap(X, l, STEP_SIZE, K, sys); % solve data driven Lyapunov equation
% Pi = dd_sim_linear.Pi;
% Kip1 = dd_sim_linear.Kip1;
% data_end_point = dd_sim_linear.data_matrix_end_point;
% data_int = dd_sim_linear.data_array_integral;
% 
% % calculate error matrix
% Delta_i = Pi - dd_sim_linear.Pi_lyap;
% norm_Delta_i = norm(Delta_i);
% fprintf(['|Delta_i| = ', '\n'])
% fprintf([num2str(norm(Delta_i)), '\n'])


%% --------------------------- Functions --------------------------
function generatePhyiscalDataPlot(title,x,t,LegendName)
% figure("Name",title)
% sgtitle(title,'Interpreter','latex')
subplot(3,1,1)
hold on
plot(t,x(1,:),':','DisplayName',LegendName(1),LineWidth=2)
plot(t,x(3,:),':','DisplayName',LegendName(3),LineWidth=2)
legend('Interpreter','latex',Location="best")
xlim([t(1),t(end)])
grid on

subplot(3,1,2)
hold on
plot(t,x(2,:),':','DisplayName',LegendName(2),LineWidth=2)
plot(t,x(4,:),':','DisplayName',LegendName(4),LineWidth=2)
legend('Interpreter','latex',Location="best")
xlim([t(1),t(end)])
grid on
end

function generateStateDataPlot(title,x,t)
figure("Name",title)
sgtitle(title,'Interpreter','latex')
plot(t,x,'.-',MarkerSize=10);
xlabel('$t$','Interpreter','latex')
ylabel('$x$','Interpreter','latex')
xlim([t(1),t(end)])
grid on
ax = gca;
ax.FontSize = 14;
end