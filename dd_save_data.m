% consatantly comment out every line
% ONLY un-comment out when use
% clear variables
% close all
% clc
%% Save Experiment Raw data - after initialization and simulink experiment
% date = "240729";
% filename = "data_raw_2407291645";
% comment = ...
%     "Usage: Data-driven data, 20 different IC. " + ...
%     "K: K1. " + ...
%     "Initial condition: case20 " + ...
%     "Differentiator: no filtering. ";
% save("data\data_" + date + "\" + filename + ".mat")

%% Save Data to Temperary file - when doing data processing
% save("data\data_temp.mat","X","T","K","STEP_SIZE")

%% Save Data-driven Data - after data processing
date = "240729";
filename = "nonlinear_sim_X_T_2407291757";
comment = "20 IC x 2 seg: [0.2,0.7] [0.7,1.2]. Raw data: all data from data_raw_240729xxxx";
save("data\data_" + date + "\" + filename + ".mat",...
    "X_sim", "T", "K", "STEP_SIZE","comment")

%% Save Simulation Raw Data
% date = "240719";
% filename = "SIM_X_T_2407192157";
% comment = "Simulation data 4 IC x 3 seg, 0.2s interval.";
% 
% save("data\data_" + date + "\" + filename + ".mat",...
%     "X", "T", "K", "STEP_SIZE","comment")