clear variables
close all
clc

%%
date = "250202";
step = 8; % CHECK BUFORE RUN!!!!!

% load all data
dir_path = "data\data_" + date + "\Step_" + num2str(step) + "\";
matlist = dir(dir_path + 'RawData*.mat');

% compute avg of K_repeat
sum_Ki_repeat = [];
sum_Pi_repeat = zeros(4,4);
for r = 1:size(matlist,1)
    load(dir_path + matlist(r).name)
    sum_Ki_repeat(r,:) = eps_Kip1;
    sum_Pi_repeat = sum_Pi_repeat + eps_Pi;
end

Kip1 = sum(sum_Ki_repeat,1)/size(matlist,1)
Pi = sum_Pi_repeat./size(matlist,1)

Delta_Ki = Kip1 - K_opt;
norm_Delta_Ki = norm(Delta_Ki);
fprintf(['|K_i+1 - K^*| = ', '\n'])
fprintf([num2str(norm_Delta_Ki), '\n'])

save(dir_path + "Solution_Step_" + num2str(step) + ".mat", "Pi","Kip1")