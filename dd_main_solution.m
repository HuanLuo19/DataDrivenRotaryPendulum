clear variables
close all
clc
%%
load("data\data_241020\X_T_2410202156.mat")
% load("data\data_240729\linear_sim_X_T_2407291757.mat")
% X_linear = X_sim;
figure("Name","All Data")
sgtitle("All Data",'Interpreter','latex')
for i = 1:size(X,1)
    subplot(size(X,1)/10,10,i) % [change layout as needed]
    % plot(T{i},X_linear{i},'r--',LineWidth=2)
    % hold on
    plot(T{i},X{i},'g-',LineWidth=1)
    hold off
    % xlabel('$t$','Interpreter','latex')
    % ylabel('$x$','Interpreter','latex')
    xlim([T{i}(1),T{i}(end)])
    grid on
end

% load("data\data_241015\linear_sim_X_T_2410152352.mat")
% X_linear = X_sim(1:40,1);
% load("data\data_241015\nonlinear_sim_X_T_2410152352.mat")
% X_nonlinear = X_sim(1:40,1);
% 
% figure("Name","All Data")
% sgtitle("All Data - Red: linear, Blue: Non-linear",'Interpreter','latex')
% for i = 1:size(X_linear,1)
%     subplot(size(X_linear,1)/10,10,i) % [change layout as needed]
%     plot(T{i},X_linear{i},'r-',LineWidth=2)
%     hold on
%     plot(T{i},X_nonlinear{i},'b--',LineWidth=2)
%     hold off
%     % xlabel('$t$','Interpreter','latex')
%     % ylabel('$x$','Interpreter','latex')
%     xlim([T{i}(1),T{i}(end)])
%     grid on
% end

%% Data for DD

% for i = 1:40
%     X(i,1) = X_linear(i,1);
% end
% for i = 1:40
%     X(i,1) = X_nonlinear(i,1);
% end
% X = [X_nonlinear(1,1); X_nonlinear(10,1); X_nonlinear(13,1); X_nonlinear(20,1); X_nonlinear(22,1); 
%     X_nonlinear(27,1); X_nonlinear(32,1); X_nonlinear(36,1); X_nonlinear(37,1); X_nonlinear(38,1)]
% X = [X_nonlinear(31:40,1)];

figure("Name","All Data")
sgtitle("All Data",'Interpreter','latex')
for i = 1:size(X,1)
    subplot(size(X,1)/10,10,i) % [change layout as needed]
    plot(T{i},X{i},'-')
    % xlabel('$t$','Interpreter','latex')
    % ylabel('$x$','Interpreter','latex')
    xlim([T{i}(1),T{i}(end)])
    grid on
end

%% Data-driven Solution
% system
A = [0    1.0000         0         0;
     0  -12.2135         0         0;
     0         0         0    1.0000;
     0   -7.6602  -66.8782   -0.2289];
B = [0   39.2743         0   24.6327]';
Q = diag([100 0 1000 0]);
R = 100;

if ~all(eig(A - B * K)<0)
    fprintf("K is NOT a stabilizing gain! \n")
    return
end 

x0 = zeros(size(A,1),1);
sys = linearSys(A,B,x0,Q,R);

%
l = 0;
% X = X{1:10,1};
dd = ddLyap(X, l, STEP_SIZE, K, sys); % solve data driven Lyapunov equation
Pi = dd.Pi;
Kip1 = dd.Kip1;
Pi_lyap = dd.Pi_lyap;
data_end_point = dd.data_matrix_end_point;
data_int = dd.data_array_integral;

Delta_i = Pi - Pi_lyap;
norm_Delta_i = norm(Delta_i);
fprintf(['|Delta_i| = ', '\n'])
fprintf([num2str(norm_Delta_i), '\n'])
if ~all(eig(A - B * Kip1)<0)
    fprintf("K_{i+1} is NOT a stabilizing gain! \n")
    fprintf("eigenvalues of A-BK_{i+1} \n ")
    eig(A - B * Kip1)
    % return
end

%% ddLyap data matrix

dd_linear = ddLyap(X_linear(1:40,1), l, STEP_SIZE, K, sys); % solve data driven Lyapunov equation
I_xx_linear = dd_linear.data_array_integral;
delta_xx_linear = dd_linear.data_matrix_end_point;
Pi_linear = dd_linear.Pi;


dd_nonlinear = ddLyap(X, l, STEP_SIZE, K, sys); % solve data driven Lyapunov equation
I_xx_nonlinear = dd_nonlinear.data_array_integral;
delta_xx_nonlinear = dd_nonlinear.data_matrix_end_point;
Pi_nonlinear = dd_nonlinear.Pi;

error_delta_xx = abs((delta_xx_linear - delta_xx_nonlinear) ./ delta_xx_linear);
error_I_xx = abs((I_xx_linear - I_xx_nonlinear) ./ I_xx_linear);
error_Pi = abs((Pi_linear - Pi_nonlinear));

% heatmap
figure("Name","Data matrix comparison")
sgtitle("Data matrix comparison: Linear Sim & Non-linear Sim - replace data No.6")
subplot(3,1,1)
heatmap(error_delta_xx)
title('Error of \delta_{xx}')
subplot(3,1,2)
heatmap(error_I_xx)
title("Error of I_{xx}")
subplot(3,1,3)
heatmap(error_Pi)
title("Error of P_i")
