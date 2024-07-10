clear variables
close all
clc

%% Initialize simulation solver
STEP_SIZE = 5e-3;
SIMULATION_TIME = inf;

Tf1 = 0.08; % derivative transfer fcn coeff
Tf2 = 0.08;
%%
K0 = [10 1 1 -0.1];
% K1 = [2086.16114717546	358.022891138783	-99.9245069315659	-20.3349766266169];

A = [0    1.0000         0         0;
     0  -12.2135         0         0;
     0         0         0    1.0000;
     0   -7.6602  -66.8782   -0.2289];
B = [0;
   39.2743;
         0;
   24.6327;];
Q = [100         0           0           0;
     0           0           0           0;
     0           0        1000           0;
     0           0           0           0];
R = 1;
disp('***** eigenvalues of (A + BK0) *****')
eig(A - B*K0)
K_opt = lqr(A,B,Q,R);
disp('***** eigenvalues of (A + BK*) *****')
eig(A - B*K_opt)

