clear variables
close all
clc

SIMULATION_TIME = 20;
STEP_SIZE = 5e-3; % step size
%% Systen Parameters
identified_para = matfile("data/identified_parameters.mat");
a1 = identified_para.a1;
a2 = identified_para.a2;
b1 = identified_para.b1;
b2 = identified_para.b2;

g  = 9.81;
L1 = 9.20e-02;  % 9.2cm

%% Initial Condition
x0 = [2 0 1 0]';
K = [10 1 10 1];

nonlinear_sim = sim('sim_nonlinear_pend.slx');

%% Plot
% t = 0:STEP_SIZE:SIMULATION_TIME;
figure("Name","Nonlinear Simulation")
sgtitle("Nonlinear Simulation",'Interpreter','latex')
% subplot(3,1,1)
hold on
plot(nonlinear_sim.tsim,nonlinear_sim.theta1sim,'-',LineWidth=2)
plot(nonlinear_sim.tsim,nonlinear_sim.phi2sim,'-',LineWidth=2)
legend("motor angle SIM","pendulum angle SIM", ...
    'Interpreter','latex',Location="best")
xlim([nonlinear_sim.tsim(1),nonlinear_sim.tsim(end)])
grid on

