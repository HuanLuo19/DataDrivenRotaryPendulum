clear variables
close all
clc

%%

date = "250202";
figure
hold on

%% system
A = [0    1.0000         0         0;
     0  -12.2135         0         0;
     0         0         0    1.0000;
     0   -7.6602  -66.8782   -0.2289];
B = [0   39.2743         0   24.6327]';
Q = diag([100 0 1000 0]);
R = 100;

x0 = zeros(size(A,1),1);
sys = linearSys(A,B,x0,Q,R);
K_opt = sys.K_opt;

P0 = [0    20   0     10;
          20   3    0     1;
          0    0    0     100;
          10   1    100   1];
K = R \ B' * P0;


% load all data
itr_Ki = K;
itr_Delta_Ki = norm(K - K_opt);
for step = 1:8
    dir_path = "data\data_" + date + "\Step_" + num2str(step) + "\";
    load(dir_path + "Solution_Step_" + num2str(step) + ".mat")
    itr_Ki = [itr_Ki; Kip1];
    itr_Delta_Ki = [itr_Delta_Ki; norm(Kip1 - K_opt)];
end



%%

for step = 1:8
    
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
        h1 = plot(step+1,norm(eps_Kip1 - K_opt),'bx',MarkerSize=12,LineWidth=1.2);
    end
    
    Kip1 = sum(sum_Ki_repeat,1)/size(matlist,1)
    Pi = sum_Pi_repeat./size(matlist,1)
    
    Delta_Ki = Kip1 - K_opt;
    norm_Delta_Ki = norm(Delta_Ki);
    fprintf(['|K_i+1 - K^*| = ', '\n'])
    fprintf([num2str(norm_Delta_Ki), '\n'])

end
h2 = plot(1:step + 1,itr_Delta_Ki,'g.-',MarkerSize=26,LineWidth=1.2);

legend([h1,h2],"$\| K_i - K^* \| $ Candidates","$\| K_i - K^* \|$ Average","Interpreter","latex")
xlabel("$i$","Interpreter","latex")
xlim([1 step + 1])
grid on
box on
ax = gca;
% ax.FontSize = 14;
ax.XTick = 1: step + 1;