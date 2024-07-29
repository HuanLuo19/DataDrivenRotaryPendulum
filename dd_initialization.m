clear variables
close all
clc
%% Initialize simulation
% set step size and simulation time
STEP_SIZE = 5e-3;
SIMULATION_TIME = inf;

% select controller gain
% K = [10 1 1 -0.1];
K = [10 1 10 1];
% K = [10 2.0595 26.7286 -1.8444]; % crane K opt
% K = [-0.7071 -0.8513 -17.4339 -2.0987]; % inverted pendulum K opt

%% system
A = [0    1.0000         0         0;
     0  -12.2135         0         0;
     0         0         0    1.0000;
     0   -7.6602  -66.8782   -0.2289];
B = [0   39.2743         0   24.6327]';
Q = diag([100 0 1000 0]);
R = 1;

if ~all(eig(A - B * K)<0)
    fprintf("K is NOT a stabilizing gain! \n")
    return
end 
Pi_lyap = lyap((A - B * K)', Q + K' * R * K);

load("data\K1_2407181636.mat")
K = Kip1;
% --- run simulink --->
%% ---> Save  Experiment Raw Data
% run "dd_save_data.m": Save Experiment Raw Data section
