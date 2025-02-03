clear variables
close all
clc

SIMULATION_TIME = 1;
STEP_SIZE = 5e-3;

%% System Parameters
identified_para = matfile("data/identified_parameters.mat");
a1 = identified_para.a1;
a2 = identified_para.a2;
b1 = identified_para.b1;
b2 = identified_para.b2;

g  = 9.81;
L1 = 9.20e-02;  % 9.2cm

A = [0    1.0000         0         0;
     0  -12.2135         0         0;
     0         0         0    1.0000;
     0   -7.6602  -66.8782   -0.2289];
B = [0   39.2743         0   24.6327]';
Q = diag([100 0 1000 0]);
R = 100;

%% Initial Condition
load("data\data_common\x0_all_2410202156.mat")

%% Initial Gain
K = [10 1 10 1];
K0 = K;

%% ----------------- ADP --------------
itr = 1;
Ki_itr = cell(1,itr);
Pi_itr = cell(1,itr);
Ki_itr{1} = K;
Ki_lyap_itr = cell(1,itr);
Pi_lyap_itr = cell(1,itr);

for j = 1:itr
    fprintf(['----------  itr = ', num2str(j), ' ------------\n'])
   
    % Data and Plot
    X = cell(0,0);         
    T = cell(0,0);            
    for i = 1:size(x0_all,2)
        x0 = x0_all(:,i);
    
        % % linear oneshot
        run("sim_linearSys_oneshot.m")
        dataPros_SimLinear2Datadriven = dataProcessing(x_sim_linear,t_sim_linear,[]);
        T_int = [0,1];
        [x_dd, t_dd, ~] = dataPros_SimLinear2Datadriven.getTimeIntervalData(T_int);
        X = [X; x_dd];
        T = [T; t_dd];
    end

    segment_number = 10;
    figure("Name","All Data-driven Data")
    sgtitle("All Data-driven Data - Linear System",'Interpreter','latex')
    for i = 1:size(X,1)
        subplot(size(X,1)/segment_number,segment_number,i)
        plot(T{i},X{i},'.-',MarkerSize=10)
        xlabel('$t$','Interpreter','latex')
        ylabel('$x$','Interpreter','latex')
        xlim([T{i}(1),T{i}(end)])
        grid on
    end

    % Data-driven Solution
    if ~all(eig(A - B * K)<0)
        fprintf("K is NOT a stabilizing gain! \n")
        return
    end 
    x0 = zeros(size(A,1),1);
    sys = linearSys(A,B,x0,Q,R);
    
    l = 0;
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
   
    e = 0.5;
    Kip1_hat = (1 - e) * K + e * Kip1;
    K = Kip1_hat
    % K = Kip1

    Ki_itr{j + 1} = K;
    Pi_itr{j} = Pi;
    Ki_lyap_itr{j + 1} = (1 - e) * K + e * dd.Kip1_lyap;
    Pi_lyap_itr{j} = dd.Pi_lyap;
    if ~all(eig(A - B * K)<0)
        fprintf("K_{i+1} is NOT a stabilizing gain! \n")
        fprintf("eigenvalues of A-BK_{i+1} \n ")
        eig(A - B * K)
        return
    end
end

%% Plot iteration result
figure
Ki_itr_plot = zeros(itr,4);
Ki_lyap_itr_plot = zeros(itr,4);
for i = 1:itr + 1
    Ki_itr_plot(i,:) = Ki_itr{i};
end
for i = 2:itr + 1
    Ki_lyap_itr_plot(i,:) = Ki_lyap_itr{i};
end
plot(0:itr, Ki_itr_plot,'.-', ...
    MarkerSize=24, LineWidth=1.5);
hold on
plot(0:itr, Ki_lyap_itr_plot,'.--', ...
    MarkerSize=24, LineWidth=1.5);
plot(itr + 1, sys.K_opt(1),'*',MarkerEdgeColor="#0072BD",MarkerSize=10, LineWidth=1.5)
plot(itr + 1, sys.K_opt(2),'*',MarkerEdgeColor="#D95319",MarkerSize=10, LineWidth=1.5)
plot(itr + 1, sys.K_opt(3),'*',MarkerEdgeColor="#EDB120",MarkerSize=10, LineWidth=1.5)
plot(itr + 1, sys.K_opt(4),'*',MarkerEdgeColor="#7E2F8E",MarkerSize=10, LineWidth=1.5)
legend('$\tilde K_i(1)$', '$\tilde K_i(2)$','$\tilde K_i(3)$','$\tilde K_i(4)$', ...
    '$K^*(1)$','$K^*(2)$','$K^*(3)$','$K^*(4)$', ...
    'Interpreter','latex',Location="best")
hold off
xlim([0,itr + 1])
xlabel('$t$','Interpreter','latex')
ylabel('Numeric value')
title('$K$ matrix elements','Interpreter','latex')
ax = gca;
ax.FontSize = 14;
%% --------------------------- Functions --------------------------
function generateStateDataPlot(x,t)
% figure("Name",title)
% sgtitle(title,'Interpreter','latex')
plot(t,x,'.-',MarkerSize=10);
% legend('$\theta_1$','$\dot \theta_1$','$\phi_2$','$\dot \phi_2$','Interpreter','latex')
% xlabel('$t$','Interpreter','latex')
% ylabel('$x$','Interpreter','latex')
xlim([t(1),t(end)])
grid on
% ax = gca;
% ax.FontSize = 14;
end
