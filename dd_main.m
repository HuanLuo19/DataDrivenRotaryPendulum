clear variables
close all
clc

%% Initialize
STEP_SIZE = 5e-3;
SIMULATION_TIME = inf;
% load("data/data_K0_2407081823.mat")
load("data/data_K0_2407101229_transFcnDeriv.mat")

%% system
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


K0 = [10 1 1 -0.1];
if ~all(eig(A - B * K0)<0)
    fprintf("K0 is NOT a stabilizing gain! \n")
    return
end 
Pi_lyap = lyap((A - B * K0)', Q + K0' * R * K0);
%% Show Raw Data
figure
sgtitle("Raw data",'Interpreter','latex')
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
legend("switch state",...
    'Interpreter','latex',Location="northwest")

%% Get data-driven data
time_on_idx = find(switch_state.Data(1,:),1); % get switch on time
l = 10; % segment(equation) numbers
T_int = 0.05;   % integral time
N_data_points = l * T_int / STEP_SIZE;

dd_motor_angle = motor_angle.Data(time_on_idx:time_on_idx + N_data_points);
dd_motor_speed = motor_speed.Data(time_on_idx:time_on_idx + N_data_points);
dd_pendulum_angle = pendulum_angle.Data(time_on_idx:time_on_idx + N_data_points);
dd_pendulum_speed = pendulum_speed.Data(time_on_idx:time_on_idx + N_data_points);
dd_time = time.Data(time_on_idx:time_on_idx + N_data_points);
dd_switch_state = switch_state.Data(1,time_on_idx:time_on_idx + N_data_points);

figure
sgtitle("Data-driven data",'Interpreter','latex')
subplot(3,1,1)
hold on
plot(dd_time,dd_motor_angle,'.-',MarkerSize=10)
plot(dd_time,dd_pendulum_angle,'.-',MarkerSize=10)
legend("motor angle","pendulum angle", ...
    'Interpreter','latex',Location="northwest")
xlim([dd_time(1),dd_time(end)])
subplot(3,1,2)
hold on
plot(dd_time,dd_motor_speed,'.-',MarkerSize=10)
plot(dd_time,dd_pendulum_speed,'.-',MarkerSize=10)
legend("motor speed","pendulum speed",...
    'Interpreter','latex',Location="northwest")
xlim([dd_time(1),dd_time(end)])
subplot(3,1,3)
hold on
plot(dd_time,dd_switch_state,'.-',MarkerSize=10)
legend("switch state",...
    'Interpreter','latex',Location="northwest")
xlim([dd_time(1),dd_time(end)])

% % RAW state data and dd solution
% x = [dd_motor_angle, dd_motor_speed, dd_pendulum_angle, dd_pendulum_speed]';
% x0 = x(:,1);
% sys = linearSys(A,B,x0,Q,R);
% if ~all(eig(A - B * K0)<0)
%     fprintf("Ki is NOT a stabilizing gain! \n")
%     % return
% end    
% dd = ddLyap(x, STEP_SIZE, l , K0, sys); % solve data driven Lyapunov equation
% Pi = dd.Pi;
% Kip1 = dd.Kip1;
% 
% Delta_i = Pi - Pi_lyap;
% fprintf(['|Delta_i| = ', '\n'])
% fprintf([num2str(norm(Delta_i)), '\n'])
% if ~all(eig(A - B * Kip1)<0)
%     fprintf("K_{i+1} is NOT a stabilizing gain! \n")
%     fprintf("eigenvalues of A-BK_{i+1} \n ")
%     eig(A - B * Kip1)
%     % return
% end   
% 
% figure
% plot(dd_time',x,'.-',MarkerSize=10);
% xlabel('$t$','Interpreter','latex')
% ylabel('$x$','Interpreter','latex')
% xlim([dd_time(1),dd_time(end)])
% ax = gca;
% ax.FontSize = 14;

% state data and dd solution
x = [dd_motor_angle, dd_motor_speed, dd_pendulum_angle, dd_pendulum_speed]';

% smoothing 
% x = smoothdata(x,2,"gaussian",20);
% smooth_mehod = "gaussian";
% winsize = 10;
% x(2,:) = smoothdata(x(2,:),num2str(smooth_mehod),winsize);
% x(4,:) = smoothdata(x(4,:),num2str(smooth_mehod),winsize);

x0 = x(:,1);
sys = linearSys(A,B,x0,Q,R);
if ~all(eig(A - B * K0)<0)
    fprintf("Ki is NOT a stabilizing gain! \n")
    % return
end    
dd = ddLyap(x, STEP_SIZE, l , K0, sys); % solve data driven Lyapunov equation
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

figure
plot(dd_time',x,'.-',MarkerSize=10);
xlabel('$t$','Interpreter','latex')
ylabel('$x$','Interpreter','latex')
xlim([dd_time(1),dd_time(end)])
ax = gca;
ax.FontSize = 14;