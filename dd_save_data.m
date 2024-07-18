% consatantly comment out every line
% ONLY un-comment out when use
% clear variables
% close all
% clc
%% Save Raw data - after initialization and simulink experiment
% date = "240718";
% filename = "data_raw_2407181645";
% comment = ...
%     "Usage: Data-driven data" + ...
%     "K: K1_2407181636. " + ...
%     "Initial condition: motor angle = 2. " + ...
%     "Differentiator: no filtering. ";
% save("data\data_" + date + "\" + filename + ".mat")

%% Save Data to Temperary file - when doing data processing
% save("data\data_temp.mat","X","T","K","STEP_SIZE")

%% Save Data-driven Data - after data processing
% date = "240718";
% filename = "X_T_2407181617";
% comment = "Raw data: data_raw_2407181539,2407181542,2407181544,2407181546";
% save("data\data_" + date + "\" + filename + ".mat",...
%     "X", "T", "K", "STEP_SIZE","comment")

