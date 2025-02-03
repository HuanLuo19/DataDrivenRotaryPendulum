%% RUN 2 TIMES TO SAVE ALL !!!!
%% RUN 2 TIMES TO SAVE ALL !!!!
%% RUN 2 TIMES TO SAVE ALL !!!!
%% Save Experiment Raw data - after initialization and simulink experiment
% date = "250202";
% filename = "RawData_250202";
% 
% comment = ...
%     "Usage: Data-driven data, 20 different IC. " + ...
%     "Initial condition: case20 " + ...
%     "Differentiator: filtering s/0.08s+1. ";
% save("data\data_" + date + "\Step_" + num2str(step) + "\"+ ...          % folder name
%     filename + "_" + num2str(step) + "." + num2str (repeat) + ".mat")   % file name
% 
% close all
% clc
%% Load data if needed
clear variables
close all
clc

date = "250202";
filename = "RawData_250202";
step = 8;
repeat = 2; % change segment_number below if [step,repeat] = [4,2] [7,4] [7,5] [8,4] [8,5]

load("data\data_" + date + "\Step_" + num2str(step) + "\" + ...
    filename + "_" + num2str(step) + "." + num2str (repeat) + ".mat")

%% Get Data-driven Dataset
% plot raw data: state, switch

legend_name_experiment = ["motor angle","motor speed","pendulum angle","pendulum speed"];
generatePhyiscalDataPlot("Raw Data",state4dim.Data',time.Data,switch_state.Data(1,:),legend_name_experiment)

% divide segments
idx_start = [];
idx_end = [];
k_start = 1;
k_end = 1;
for i = 1:size(time.data,1) - 1 
    if switch_state.Data(1,i) == 0 && switch_state.Data(1,i+1) == 1
        idx_start(k_start) = i + 1;
        k_start = k_start + 1;
    end
    if switch_state.Data(1,i) == 1 && switch_state.Data(1,i+1) == 0
        idx_end(k_end) = i;
        k_end = k_end + 1;
    end
end
segment_number = size(idx_start,2);
% segment_number = 20; % for step 4, repeat 2
% segment_number = 19; % for step 7, repeat 4
% segment_number = 13; % for step 7, repeat 5
% segment_number = 14; % for step 8, repeat 4
% segment_number = 19; % for step 8, repeat 5

% construct raw data set
X_raw = cell(segment_number,1);
T_raw = cell(segment_number,1);
for i = 1:segment_number
    x_raw = state4dim.Data(idx_start(i):idx_end(i),:)';
    t_raw = time.Data(idx_start(i):idx_end(i)) - time.Data(idx_start(i));
    X_raw{i} = x_raw;
    T_raw{i} = t_raw;
end
figure("Name","State segments - Raw")
sgtitle("State segments - Raw",'Interpreter','latex')
for i = 1:segment_number
    subplot(4,5,i)
    plot(T_raw{i},X_raw{i},'.-',MarkerSize=10)
    xlabel('$t$','Interpreter','latex')
    ylabel('$x$','Interpreter','latex')
    xlim([T_raw{i}(1),T_raw{i}(end)])
    grid on
end

% get dd data: set time interval from raw data set
T_int = [0.5,2.5];
X = cell(segment_number,1);
T = cell(segment_number,1);
for i = 1:segment_number
    dd_x = X_raw{i}(:,T_int(1)/STEP_SIZE:T_int(2)/STEP_SIZE);
    dd_t = T_raw{i}(T_int(1)/STEP_SIZE:T_int(2)/STEP_SIZE);
    X{i} = dd_x;
    T{i} = dd_t;
end
figure("Name","State segments")
sgtitle("State segments",'Interpreter','latex')
for i = 1:segment_number
    subplot(4,5,i)
    plot(T{i},X{i},'.-',MarkerSize=10)
    xlabel('$t$','Interpreter','latex')
    ylabel('$x$','Interpreter','latex')
    xlim([T{i}(1),T{i}(end)])
    grid on
end

%% Data-driven Solution
% system
A = [0    1.0000         0         0;
     0  -12.2135         0         0;
     0         0         0    1.0000;
     0   -7.6602  -66.8782   -0.2289];
% B = [0   39.2743         0   24.6327]';
% Q = diag([100 0 1000 0]);
% R = 100;

if ~all(eig(A - B * K)<0)
    fprintf("K is NOT a stabilizing gain! \n")
    return
end 

x0 = zeros(size(A,1),1);
sys = linearSys(A,B,x0,Q,R);

% ddlyap solution
l = 0;
epsilon = 0.01;
eps_dd = eps_ddLyap(X, l, STEP_SIZE, P0, K, sys, epsilon);

% ---Pi---
eps_Pi = eps_dd.Pi_e;
eps_Pi_lyap = eps_dd.Pi_lyap_e;
eps_Delta_Pi = eps_Pi - eps_Pi_lyap;
eps_norm_Delta_Pi = norm(eps_Delta_Pi);
% fprintf(['|eps_Delta_Pi| = ', '\n'])
% fprintf([num2str(eps_norm_Delta_Pi), '\n'])

% ---Ki---
K
eps_Kip1 = eps_dd.Kip1_e
eps_Kip1_lyap = eps_dd.Kip1_lyap_e
K_opt = sys.K_opt;
eps_Delta_K = K - K_opt;
eps_norm_Delta_K = norm(eps_Delta_K);
fprintf(['|K_i - K^*| = ', '\n'])
fprintf([num2str(eps_norm_Delta_K), '\n'])

eps_Delta_Ki = eps_Kip1 - K_opt;
eps_norm_Delta_Ki = norm(eps_Delta_Ki);
fprintf(['|K_i+1,r - K^*| = ', '\n'])
fprintf([num2str(eps_norm_Delta_Ki), '\n'])

if ~all(eig(A - B * eps_Kip1)<0)
    fprintf("K_{i+1} is NOT a stabilizing gain! \n")
    fprintf("eigenvalues of A-BK_{i+1} \n ")
    eig(A - B * eps_Kip1)
    % return
end

% % for calibration
% l = 0;
% dd = ddLyap(X, l, STEP_SIZE, K, sys); % solve data driven Lyapunov equation
% Pi = dd.Pi;
% Kip1 = dd.Kip1;
% Pi_lyap = dd.Pi_lyap;
% 
% data_end_point = dd.data_matrix_end_point;
% data_int = dd.data_array_integral;
% 
% Delta_Pi = Pi - Pi_lyap;
% norm_Delta_Pi = norm(Delta_Pi);
% fprintf(['|Delta_Pi| = ', '\n'])
% fprintf([num2str(norm_Delta_Pi), '\n'])
% if ~all(eig(A - B * Kip1)<0)
%     fprintf("K_{i+1} is NOT a stabilizing gain! \n")
%     fprintf("eigenvalues of A-BK_{i+1} \n ")
%     eig(A - B * Kip1)
%     % return
% end
% 
% eps_Pi_lyap - (P0 * (1- epsilon) + epsilon * Pi_lyap)

%% Generate Plots Fcns
function generatePhyiscalDataPlot(title,x,t,swt,LegendName)
    figure("Name",title)
    sgtitle(title,'Interpreter','latex')
    subplot(3,1,1)
    hold on
    plot(t,x(1,:),'.-','DisplayName',LegendName(1),MarkerSize=10)
    plot(t,x(3,:),'.-','DisplayName',LegendName(3),MarkerSize=10)
    legend('Interpreter','latex',Location="best")
    xlabel('$t$','Interpreter','latex')
    xlim([t(1),t(end)])
    grid on
    
    subplot(3,1,2)
    hold on
    plot(t,x(2,:),'.-','DisplayName',LegendName(2),MarkerSize=10)
    plot(t,x(4,:),'.-','DisplayName',LegendName(4),MarkerSize=10)
    legend('Interpreter','latex',Location="best")
    xlabel('$t$','Interpreter','latex')
    xlim([t(1),t(end)])
    grid on
    
    subplot(3,1,3)
    hold on
    plot(t,swt,'.-',MarkerSize=10)
    legend("switch state",...
        'Interpreter','latex',Location="best")
    xlabel('$t$','Interpreter','latex')
    xlim([t(1),t(end)])
end