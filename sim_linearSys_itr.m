clear variables
close all
clc

%% Define systen
identified_para = matfile("data/identified_parameters.mat");
a1 = identified_para.a1;
a2 = identified_para.a2;
b1 = identified_para.b1;
b2 = identified_para.b2;

g  = 9.81;
L1 = 9.20e-02;  % 9.2cm

E = [1 0 0 0;
     0 1 0 0;
     0 0 1 0;
     0 -L1 0 a2 ];
F = [ 0  1  0  0
      0 -a1 0  0
      0  0  0  1
      0  0 -g -b2 ];
G = [ 0
      b1
      0
      0 ];
A = E \ F;
B = E \ G;
% B = [1 1 1 1]';
nx = size(A,2); %state dimension
nu = size(B,2); %input dimension
q1 = 100;
q3 = 1000;
Q = diag([q1 0 q3 0]);
R = 1;

x0 = [1 0 0 0]';
sys_0 = linearSys(A,B,x0,Q,R);

%% TEST
% test1
% T_test = 5;
% sys_0.generateTestResponseZeroInput(T_test);
% sys_0.generateTestResponseOptimalInput(T_test);
% 
% test2
% dtau = 1e-3;
% T = 0.5;
% K_a = -[-1 -0.1 -0.1 0.1];
% [x_a,u_a,tspan_a] = sys_0.ResponseFromGain(K_a,dtau,T);
% [x_opt,u_opt,tspan_opt] = sys_0.ResponseFromGain(K_opt,dtau,T);

%% Data loss processing
% T_loss = 0.05;
% i_loss = 1;
% loss1 = dataLoss(x_a,tspan_a,i_loss,T_loss);
% x_a_loss = loss1.x_loss;
% loss1.showLossData()
% loss1.generateLinearInterpolatedData()
% x_a_intrp = loss1.x_intrp;

%% ADP
dtau = 5e-3;
l = 10;         % number of equations(data segments)
T_int = 0.05;   % integral time
T = l * T_int;  % iteration time

% K0 = [10 1 1 -0.1];
K0 = [10 1 10 1];
% K0 = [330.7524   34.1486  435.6918  -33.4793];
% K1 = [2086.16114717546	358.022891138783	-99.9245069315659	-20.3349766266169];
% K1 = [1693.85882885054	317.304940008585	-471.780137355098	-68.7236527650644];
% K0 = K1;

if ~all(eig(A - B * K0)<0)
    fprintf("K0 is NOT a stabilizing gain! \n")
    return
end
itr = 30;

Flag_Sim_No_Loss = 1;
Flag_Sim_Loss = 0;

%% Sim Loss -----------------------------------------
% if Flag_Sim_Loss
% % loss data
% x0_loss = x0;
% K0_loss = K0;
% x_itr_loss = [];
% u_itr_loss = [];
% t_itr_loss = [];
% Ki_itr_loss = cell(1,itr);
% Pi_itr_loss = cell(1,itr);
% Ki_itr_loss{1} = K0;
% 
% % error data
% Delta_i_itr = cell(1,itr);
% norm_Delta_i_itr = zeros(1,itr);
% lambda_max_AmBK_itr = zeros(1,itr + 1);
% lambda_max_AmBK_itr(1) = [max(real(eig(A - B*K0_loss)))];
% for i = 1:itr
%     fprintf(['----------  itr = ', num2str(i), ' ------------\n'])
% 
%     % -------------------- loss --------------------------
%     sys_loss = linearSys(A,B,x0_loss,Q,R);
%     [x_loss_raw,u_loss,tspan] = sys_loss.ResponseFromGain(K0_loss,dtau,T);
%     if ~all(eig(A - B * K0_loss)<0)
%         fprintf("Ki_loss is NOT a stabilizing gain! \n")
%         return
%     end
%     % ----------- loss parameter
%     % T_loss = 0.5;
%     i_loss = 1;
%     loss = dataLoss(x_loss_raw,tspan,i_loss,T_loss);
%     x_loss = loss.x_loss;
%     loss.showLossData()
%     loss.generateLinearInterpolatedData()
%     x_intrp = loss.x_intrp;
%     % ----------- 
%     dd_loss = ddLyap(x_intrp, l, dtau , K0_loss, sys_loss); % solve data driven Lyapunov equation
%     Pi_loss = dd_loss.Pi;
%     Kip1_loss = dd_loss.Kip1;
%     K0_loss = dd_loss.Kip1;
%     x0_loss = x_intrp(:,end);
%     % log data
%     x_itr_loss = [x_itr_loss x_intrp];
%     u_itr_loss = [u_itr_loss u_loss];
%     t_itr_loss = [t_itr_loss tspan + T*(i-1)];
%     Ki_itr_loss{i + 1} = dd_loss.Kip1;
%     Pi_itr_loss{i} = dd_loss.Pi;
% 
%     % calculate error matrix
%     Delta_i = Pi_loss - dd_loss.Pi_lyap;
%     Delta_i_itr{i} = Delta_i;
%     norm_Delta_i_itr(i) = norm(Delta_i);
%     fprintf(['|Delta_i| = ', '\n'])
%     fprintf([num2str(norm(Delta_i)), '\n'])
% 
%     lambda_max_AmBK_itr(i + 1) = [max(real(eig(A - B*Kip1_loss)))];
% end
% % % plot
% % figure
% % subplot(2,1,1)
% % plot(t_itr_loss,x_itr_loss,'.',MarkerSize=10);
% % xlabel('$t$','Interpreter','latex')
% % title('$x-loss$','Interpreter','latex')
% % ax = gca;
% % ax.FontSize = 14;
% % 
% % subplot(2,1,2)
% % plot(t_itr_loss,u_itr_loss,'.',MarkerSize=10);
% % xlabel('$t$','Interpreter','latex')
% % title('$u$','Interpreter','latex')
% % ax = gca;
% % ax.FontSize = 14;
% end

%% Sim No Loss-----------------------------------
if Flag_Sim_No_Loss
% no loss data
x_itr = [];
u_itr = [];
t_itr = [];
Ki_itr = cell(1,itr);
Pi_itr = cell(1,itr);
Ki_itr{1} = K0;

% error data
Delta_i_itr = cell(1,itr);
norm_Delta_i_itr = zeros(1,itr);
lambda_max_AmBK_itr = zeros(1,itr);
lambda_max_AmBK_itr(1) = [max(real(eig(A - B*K0)))];
for i = 1:itr
    fprintf(['----------  itr = ', num2str(i), ' ------------\n'])
    
    % ------------------- no loss --------------------
    sys = linearSys(A,B,x0,Q,R);
    [x,u,tspan] = sys.ResponseFromGain(K0,dtau,T);
    if ~all(eig(A - B * K0)<0)
        fprintf("Ki is NOT a stabilizing gain! \n")
        return
    end
    
    % dd = ddLyap(x, l, dtau, K0, sys); % solve data driven Lyapunov equation
    % --- test for seperate data dd solution ---
    X = cell(10,1);
    X{1,1} = x(:,1:11);
    X{2,1} = x(:,11:21);
    X{3,1} = x(:,21:31);
    X{4,1} = x(:,31:41);
    X{5,1} = x(:,41:51);
    X{6,1} = x(:,51:61);
    X{7,1} = x(:,61:71);
    X{8,1} = x(:,71:81);
    X{9,1} = x(:,81:91);
    X{10,1} = x(:,91:101);
    l = 0;
    dd = ddLyap(X, l, dtau, K0, sys); % solve data driven Lyapunov equation
    % --- end test for seperate data dd solution ---
    Pi = dd.Pi;
    Kip1 = dd.Kip1;

    e = 0.01;
    Kip1_hat = (1 - e) * K0 + e * Kip1;
    K0 = Kip1_hat;
    % K0 = Kip1;

    x0 = x(:,end);
    % log data
    x_itr = [x_itr x];
    u_itr = [u_itr u];
    t_itr = [t_itr tspan + T*(i-1)];
    Ki_itr{i + 1} = dd.Kip1;
    Pi_itr{i} = dd.Pi;

    % calculate error matrix
    Delta_i = Pi - dd.Pi_lyap;
    Delta_i_itr{i} = Delta_i;
    norm_Delta_i_itr(i) = norm(Delta_i);
    fprintf(['|Delta_i| = ', '\n'])
    fprintf([num2str(norm(Delta_i)), '\n'])

    lambda_max_AmBK_itr(i + 1) = [max(real(eig(A - B*Kip1)))];
end
% plot
figure
subplot(2,1,1)
plot(t_itr,x_itr,'.',MarkerSize=10);
xlabel('$t$','Interpreter','latex')
title('$x$','Interpreter','latex')
ax = gca;
ax.FontSize = 14;

subplot(2,1,2)
plot(t_itr,u_itr,'.',MarkerSize=10);
xlabel('$t$','Interpreter','latex')
title('$u$','Interpreter','latex')
ax = gca;
ax.FontSize = 14;

figure
Ki_itr_plot = zeros(itr,nx);
for i = 1:itr + 1
    Ki_itr_plot(i,:) = Ki_itr{i};
end
% %
% Ki_itr_plot = [10 1 10 1;
%     5.9934 0.5486 6.8171 0.3307;
%     4.5409 0.2752 5.4382 0.0631;
%     3.9809 0.1691 4.9928 -0.0086;
%     3.8471 0.1344 4.9866 0.0033;
%     3.9456 0.1222 5.0904 0.0047];
% %
plot(0:itr, Ki_itr_plot,'.-', ...
    MarkerSize=24, LineWidth=1.5);
hold on
plot(itr + 1, sys_0.K_opt(1),'*',MarkerEdgeColor="#0072BD",MarkerSize=10, LineWidth=1.5)
plot(itr + 1, sys_0.K_opt(2),'*',MarkerEdgeColor="#D95319",MarkerSize=10, LineWidth=1.5)
plot(itr + 1, sys_0.K_opt(3),'*',MarkerEdgeColor="#EDB120",MarkerSize=10, LineWidth=1.5)
plot(itr + 1, sys_0.K_opt(4),'*',MarkerEdgeColor="#7E2F8E",MarkerSize=10, LineWidth=1.5)
legend('$\hat K$(1)', '$\hat K$(2)','$\hat K$(3)','$\hat K$(4)', ...
    '$K^*$(1)','$K^*$(2)','$K^*$(3)','$K^*$(4)', ...
    'Interpreter','latex',Location="best")
hold off
xlim([0,itr + 1])
xlabel('$t$','Interpreter','latex')
ylabel('Numeric value')
title('$K$ matrix elements','Interpreter','latex')
ax = gca;
ax.FontSize = 14;

end