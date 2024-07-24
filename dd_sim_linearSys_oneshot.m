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

% Generate Raw Simulation Data
STEP_SIZE = 5e-3; % step size
% IC
x0 = [-2 0 0 0]';
% FB Gain
K = [10 1 10 1];
% K = [10 1 1 -0.1];
% K1 = [330.7524   34.1486  435.6918  -33.4793];
% K = K1;
if ~all(eig(A - B * K)<0)
    fprintf("K0 is NOT a stabilizing gain! \n")
    return
end
sys = linearSys(A,B,x0,Q,R);

% Generate sim data
T_sim = 10;
[x_sim_raw,u_sim_raw,t_sim_raw] = sys.ResponseFromGain(K,STEP_SIZE,T_sim);
generateStateDataPlot("Simulation State Data",x_sim_raw,t_sim_raw)
% hold on
% generatePhyiscalDataPlot("a",x,tspan) % use this to compare with experiment data

%% Select Data from Sim Raw Data
X = cell(0,0);            % comment out this line to add data
T = cell(0,0);            % comment out this line to add data
% load("data\data_temp.mat")  % comment out this line when add the first data group

% select time interval and segment number and segment offset size
dataPros_SimRaw2Datadriven = dataProcessing(x_sim_raw,t_sim_raw,[]);
T_int = [0,0.2];
offset_size = 0.05;
segment_number = 10;
for i = 1:segment_number
    [x_dd, t_dd, ~] = dataPros_SimRaw2Datadriven.getTimeIntervalData(T_int);
    X = [X; x_dd];
    T = [T; t_dd];
    T_int = T_int + offset_size;
end

%% Plot ALL Data
figure("Name","All Data")
sgtitle("All Data",'Interpreter','latex')
for i = 1:size(X,1)
    subplot(size(X,1)/segment_number,segment_number,i)
    plot(T{i},X{i},'.-',MarkerSize=10)
    xlabel('$t$','Interpreter','latex')
    ylabel('$x$','Interpreter','latex')
    xlim([T{i}(1),T{i}(end)])
    grid on
end
%% Data-driven Solution
l = 0;
dd_sim = ddLyap(X, l, STEP_SIZE, K, sys); % solve data driven Lyapunov equation
Pi = dd_sim.Pi;
Kip1 = dd_sim.Kip1;
data_end_point = dd_sim.data_matrix_end_point;
data_int = dd_sim.data_array_integral;

% calculate error matrix
Delta_i = Pi - dd_sim.Pi_lyap;
norm_Delta_i = norm(Delta_i);
fprintf(['|Delta_i| = ', '\n'])
fprintf([num2str(norm(Delta_i)), '\n'])

%%
function generatePhyiscalDataPlot(title,x,t)
% figure("Name",title)
% sgtitle(title,'Interpreter','latex')
subplot(3,1,1)
hold on
plot(t,x(1,:),'--',LineWidth=2)
plot(t,x(3,:),'--',LineWidth=2)
legend("motor angle","pendulum angle", ...
    "motor angle SIM","pendulum angle SIM", ...
    'Interpreter','latex',Location="best")
xlim([t(1),t(end)])
grid on

subplot(3,1,2)
hold on
plot(t,x(2,:),'--',LineWidth=2)
plot(t,x(4,:),'--',LineWidth=2)
legend("motor speed","pendulum speed",...
    "motor speed SIM","pendulum speed SIM",...
    'Interpreter','latex',Location="best")
xlim([t(1),t(end)])
grid on
end

function generateStateDataPlot(title,x,t)
figure("Name",title)
sgtitle(title,'Interpreter','latex')
plot(t,x,'.-',MarkerSize=10);
xlabel('$t$','Interpreter','latex')
ylabel('$x$','Interpreter','latex')
xlim([t(1),t(end)])
grid on
ax = gca;
ax.FontSize = 14;
end