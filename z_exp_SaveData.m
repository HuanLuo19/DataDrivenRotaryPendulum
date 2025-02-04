% consatantly comment out every line
% ONLY un-comment out when use
% clear variables
% close all
% clc
%% Save Experiment Raw data - after initialization and simulink experiment
% date = "241020";
% filename = "data_raw_2410201900";
% comment = ...
%     "Usage: Data-driven data, 20 different IC. " + ...
%     "K: K1 = [10 1 10 1]. " + ...
%     "Initial condition: case20 " + ...
%     "Differentiator: filtering s/0.08s+1. ";
% save("data\data_" + date + "\" + filename + ".mat")

%% Save Data to Temperary file - when doing data processing
% save("data\data_temp.mat","X","T","K","STEP_SIZE")

%% Save Data-driven Data - after data processing
% date = "241020";
% filename = "X_T_2410211710";
% comment = "20 IC x 1 seg: [0.5,0.6]. Raw data: all data from data_raw_2410201900";
% save("data\data_" + date + "\" + filename + ".mat",...
%     "X", "T", "K", "STEP_SIZE","comment")

%% Save Simulation Data
% date = "241020";
% filename = "SIM_X_T_2410202156";
% comment = "Simulation data 20 IC x 1 seg, 1s interval. K=K1 case R=100 linear itaration";
% 
% save("data\data_" + date + "\" + filename + ".mat",...
%     "X_sim", "T", "K", "STEP_SIZE","comment")

%% Save x0
% date = "241020";
% filename = "x0_all_2410202156";
% comment = "x0 of 20 IC x 1 seg: [0.5,1.5] X_T_2410202156. Raw data: all data from data_raw_2410201900";
% save("data\data_" + date + "\" + filename + ".mat",...
%     "x0_all","comment")
