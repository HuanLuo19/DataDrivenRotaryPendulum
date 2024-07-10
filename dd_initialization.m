clear variables
close all
clc

%% Initialize simulation
% set step size and simulation time
STEP_SIZE = 5e-3;
SIMULATION_TIME = inf;

% select controller gain
% K = [10 1 1 -0.1];
% K = [10 2.0595 26.7286 -1.8444]; % crane K opt
K = [-0.7071 -0.8513 -17.4339 -2.0987]; % inverted pendulum K opt

% --- run simulink --->
%% Save data
% save("data/data_K0_2407101229_transFcnDeriv.mat")

