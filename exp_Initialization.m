clear variables
close all
clc
%% Initialize Experiment
% set step size and simulation time
STEP_SIZE = 5e-3;
SIMULATION_TIME = inf;

% Master's Thesis Setting ---
B = [0   39.2743         0   24.6327]';
Q = diag([100 0 1000 0]);
R = 100;

step = 8;
repeat = 5;

if step == 1
    P0 = [0    20   0     10;
          20   3    0     1;
          0    0    0     100;
          10   1    100   1];
    K = R \ B' * P0; % K1 = [10.3181    1.4246   24.6327    0.6391]
else
    load("data\data_250202\Step_" + num2str(step - 1) + ...
        "\Solution_Step_" + num2str(step - 1) + ".mat")
    P0 = Pi;
    K = Kip1;
end
fprintf(['K = ', '\n'])
fprintf([num2str(K), '\n'])

%% Run Experiment Simulink ---------------------------
run("DataDrivenRotaryPendulum.slx")

%%
% select controller gain
% K = [10 1 1 -0.1];
% K = [10 1 10 1]; % crane stable K
% K = [10 2.0595 26.7286 -1.8444]; % crane opt K
% K = [-0.7071 -0.8513 -17.4339 -2.0987]; % inverted opt K
% K = [-0.7071 -0.8513 -17.4339 -2.0987]; % inverted opt K

%% system
% A = [0    1.0000         0         0;
%      0  -12.2135         0         0;
%      0         0         0    1.0000;
%      0   -7.6602  -66.8782   -0.2289];
% B = [0   39.2743         0   24.6327]';
% Q = diag([100 0 1000 0]);
% R = 100;
% 
% if ~all(eig(A - B * K)<0)
%     fprintf("K is NOT a stabilizing gain! \n")
%     return
% end 
% Pi_lyap = lyap((A - B * K)', Q + K' * R * K);
% 
% load("data\K1_2407181636.mat")
% K = Kip1;
% --- run simulink --->
%% ---> Save  Experiment Raw Data
% run "exp_SaveData.m": Save Experiment Raw Data section
