clear variables
close all
clc
%%
load("test scripts\dd_linear.mat")
% dd_linear = dd;
% X_linear = X(1:16);
X_linear = X;
load("test scripts\dd_nonlinear.mat")
% dd_nonlinear = dd;
% X_nonlinear = X(1:16);
X_nonlinear = X;

A = [0    1.0000         0         0;
     0  -12.2135         0         0;
     0         0         0    1.0000;
     0   -7.6602  -66.8782   -0.2289];
B = [0   39.2743         0   24.6327]';
Q = diag([100 0 1000 0]);
R = 100;
K = [10 1 10 1];

if ~all(eig(A - B * K)<0)
    fprintf("K is NOT a stabilizing gain! \n")
    return
end 

x0 = zeros(size(A,1),1);
sys = linearSys(A,B,x0,Q,R);

%
STEP_SIZE = 5e-3;
l = 0;
dd_linear = ddLyap(X_linear, l, STEP_SIZE, K, sys);
dd_nonlinear = ddLyap(X_nonlinear, l, STEP_SIZE, K, sys);
%%
T = 0:STEP_SIZE:1;
figure("Name","All Data")
sgtitle("All Data - Red: linear, Blue: Non-linear, Green: Difference",'Interpreter','latex')
for i = 1:size(X_linear,1)
    subplot(2,10,i) % [change layout as needed]
    plot(T,X_linear{i},'r-',LineWidth=2)
    hold on
    plot(T,X_nonlinear{i},'b--',LineWidth=2)
    plot(T,X_linear{i} - X_nonlinear{i},'g-.',LineWidth=2)
    hold off
    % xlabel('$t$','Interpreter','latex')
    % ylabel('$x$','Interpreter','latex')
    xlim([T(1),T(end)])
    grid on
end

%% 验算
A1 = dd_linear.data_matrix_end_point;
b1 = - dd_linear.data_array_integral * reshape(dd_linear.Q + dd_linear.Ki' * dd_linear.R * dd_linear.Ki,[4*4,1]);
x1_lsqmn = lsqminnorm(A1,b1);
x1_pinv = pinv(A1) * b1;
x1_bkslsh = A1 \ b1;
x1 = x1_bkslsh;
eqerr1 = A1 * x1 - b1;

A2 = dd_nonlinear.data_matrix_end_point;
b2 = - dd_nonlinear.data_array_integral * reshape(dd_nonlinear.Q + dd_nonlinear.Ki' * dd_nonlinear.R * dd_nonlinear.Ki,[4*4,1]);
x2_lsqmn = lsqminnorm(A2,b2);
x2_pinv = pinv(A2) * b2;
x2_bkslsh = A2 \ b2;
x2 = x2_bkslsh;
eqerr2 = A2 * x2 - b2;

% plot
figure
subplot(2,3,1)
imagesc(A1 - A2)
colorbar
title('$A_1 - A_2$','Interpreter','latex')
subplot(2,3,2)
imagesc(b1 - b2)
title('$b_1 - b_2$','Interpreter','latex')
colorbar
subplot(2,3,3)
imagesc(x1 - x2)
title('$x_1 - x_2$','Interpreter','latex')
colorbar
subplot(2,3,4)
imagesc(eqerr1)
title('$A_1x_1 - b_1$','Interpreter','latex')
colorbar
subplot(2,3,5)
imagesc(eqerr2)
title('$A_2x_2 - b_2$','Interpreter','latex')
colorbar

%% 更新率
% P1_lyap = lyap((A - B*K)', K'*R*K + Q)
% P1 = P1_lyap;
% Ric_P1 = A' * P1 + P1 * A + Q - P1 * B * inv(R) * B' * P1;
% 
% K2 = R \ B' * P1;
% P21 = lyap((A - B*K2)', K2'*R*K2 + Q)
% P22 = P1 + lyap((A - B * inv(R) * B' * P1)', Ric_P1)
% P21 - P22
