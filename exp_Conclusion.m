clear variables
close all
clc
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

% if ~all(eig(A - B * K)<0)
%     fprintf("K is NOT a stabilizing gain! \n")
%     return
% end 
%%
date = "250202";

% load all data
itr_Ki = K;
itr_Delta_Ki = norm(K - K_opt);
for step = 1:8
    dir_path = "data\data_" + date + "\Step_" + num2str(step) + "\";
    load(dir_path + "Solution_Step_" + num2str(step) + ".mat")
    itr_Ki = [itr_Ki; Kip1];
    itr_Delta_Ki = [itr_Delta_Ki; norm(Kip1 - K_opt)];
end

% Delta_Ki = Kip1 - K_opt;
% norm_Delta_Ki = norm(Delta_Ki);
% fprintf(['|K_i+1 - K^*| = ', '\n'])
% fprintf([num2str(norm_Delta_Ki), '\n'])

plot(1:step + 1,itr_Delta_Ki,'.-',MarkerSize=20)
xlabel("$i$","Interpreter","latex")
legend("$\|K_i - K^*\|$","Interpreter","latex")
xlim([1 step + 1])
grid on
ax = gca;
ax.FontSize = 14;
ax.XTick = 1: step + 1;