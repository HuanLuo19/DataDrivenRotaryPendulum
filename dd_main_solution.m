clear variables
close all
clc
load("data\X_T_2407171746.mat")

%% Plot ALL Data
figure("Name","All Data")
sgtitle("All Data",'Interpreter','latex')
for i = 1:size(X,1)
    subplot(size(X,1)/1,1,i) % change layout as needed
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
R = 1;

if ~all(eig(A - B * K)<0)
    fprintf("K is NOT a stabilizing gain! \n")
    return
end 

x0 = zeros(size(A,1),1);
sys = linearSys(A,B,x0,Q,R);

%%
l = 0;
dd = ddLyap(X, l, STEP_SIZE, K, sys); % solve data driven Lyapunov equation
Pi = dd.Pi;
Kip1 = dd.Kip1;
Pi_lyap = dd.Pi_lyap;

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
