
%% Load Data 
% --- comment out this block 
% --- when using the date directly from the workspace
clear variables
close all
clc
% --- select data
load("data/data_K0_2407161545.mat")

%% Show Raw Data
figure("Name","Raw Data")
sgtitle("Raw Data",'Interpreter','latex')
subplot(3,1,1)
hold on
plot(motor_angle,'.-',MarkerSize=10)
plot(pendulum_angle,'.-',MarkerSize=10)
legend("motor angle","pendulum angle", ...
    'Interpreter','latex',Location="northwest")
subplot(3,1,2)
hold on
plot(motor_speed,'.-',MarkerSize=10)
plot(pendulum_speed,'.-',MarkerSize=10)
legend("motor speed","pendulum speed",...
    'Interpreter','latex',Location="northwest")
subplot(3,1,3)
hold on
plot(switch_state(1,:),'.-',MarkerSize=10)
xlabel("$t$",'Interpreter','latex')
legend("switch state",...
    'Interpreter','latex',Location="northwest")

%% Get data-driven state data
time_on_idx = find(switch_state.Data(1,:),1); % get switch on time
% l = 10; % segment(equation) numbers
% T_int = 0.05;   % integral time
% N_data_points = l * T_int / STEP_SIZE;
T_data = 1.2;
N_data_points = T_data / STEP_SIZE;

dd_time = time.Data(time_on_idx:time_on_idx + N_data_points);
dd_switch_state = switch_state.Data(1,time_on_idx:time_on_idx + N_data_points);
% for new data
dd_x = state4dim.Data(time_on_idx:time_on_idx + N_data_points,:)';
% for old data (before 240711)
% dd_x(1,:) = motor_angle.Data(time_on_idx:time_on_idx + N_data_points);
% dd_x(2,:) = motor_speed.Data(time_on_idx:time_on_idx + N_data_points);
% dd_x(3,:) = pendulum_angle.Data(time_on_idx:time_on_idx + N_data_points);
% dd_x(4,:) = pendulum_speed.Data(time_on_idx:time_on_idx + N_data_points);

% show data to use
figure("Name","Used Data")
sgtitle("Used Data",'Interpreter','latex')
subplot(3,1,1)
hold on
plot(dd_time,dd_x(1,:),'.-',MarkerSize=10)
plot(dd_time,dd_x(3,:),'.-',MarkerSize=10)
legend("motor angle","pendulum angle", ...
    'Interpreter','latex',Location="northwest")
xlim([dd_time(1),dd_time(end)])
subplot(3,1,2)
hold on
plot(dd_time,dd_x(2,:),'.-',MarkerSize=10)
plot(dd_time,dd_x(4,:),'.-',MarkerSize=10)
legend("motor speed","pendulum speed",...
    'Interpreter','latex',Location="northwest")
xlim([dd_time(1),dd_time(end)])
subplot(3,1,3)
hold on
plot(dd_time,dd_switch_state,'.-',MarkerSize=10)
legend("switch state",...
    'Interpreter','latex',Location="northwest")
xlim([dd_time(1),dd_time(end)])

figure("Name","State Data")
title('State Data','Interpreter','latex')
plot(dd_time',dd_x,'.-',MarkerSize=10);
xlabel('$t$','Interpreter','latex')
ylabel('$x$','Interpreter','latex')
xlim([dd_time(1),dd_time(end)])
ax = gca;
ax.FontSize = 14;

% smoothing 
% dd_x = smoothdata(dd_x,2,"gaussian",20);
% smooth_mehod = "gaussian";
% winsize = 10;
% dd_x(2,:) = smoothdata(dd_x(2,:),num2str(smooth_mehod),winsize);
% dd_x(4,:) = smoothdata(dd_x(4,:),num2str(smooth_mehod),winsize);

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

x0 = dd_x(:,1);
sys = linearSys(A,B,x0,Q,R);
if ~all(eig(A - B * K)<0)
    fprintf("Ki is NOT a stabilizing gain! \n")
    % return
end

l = 10;
dd = ddLyap(dd_x, l, STEP_SIZE , K, sys); % solve data driven Lyapunov equation
Pi = dd.Pi;
Kip1 = dd.Kip1;

Delta_i = Pi - Pi_lyap;
fprintf(['|Delta_i| = ', '\n'])
fprintf([num2str(norm(Delta_i)), '\n'])
if ~all(eig(A - B * Kip1)<0)
    fprintf("K_{i+1} is NOT a stabilizing gain! \n")
    fprintf("eigenvalues of A-BK_{i+1} \n ")
    eig(A - B * Kip1)
    % return
end   
