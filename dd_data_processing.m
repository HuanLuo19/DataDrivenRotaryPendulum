%% Load Data 
% --- comment out this block 
% --- when using the date directly from the workspace
clear variables
close all
clc
% --- load raw data
load("data/data_K0_2407161541.mat")
% --- Show Raw Data
generatePhyiscalDataPlot("Raw Data",state4dim.Data',time.Data,switch_state.Data(1,:))

%% Get Switch ON Data
idx_switchON_all = find(switch_state.Data(1,:)); % get all index when the switch is on

on_time = time.Data(idx_switchON_all)';
on_time = on_time - on_time(1); % RESET time label from 0
on_switch_state = switch_state.Data(1,idx_switchON_all);
on_x = state4dim.Data(idx_switchON_all,:)';

generatePhyiscalDataPlot("Switch ON Data",on_x,on_time,on_switch_state)

%% Select and Save Data from Switch ON Data
dataPros = dataProcessing(on_x,on_time);
% select time interval
T_int1 = [0,1];
T_int2 = [0.5,1.5];
T_int3 = [1,2];
[dd_x1, dd_t1] = dataPros.getTimeIntervalData(T_int1);
[dd_x2, dd_t2] = dataPros.getTimeIntervalData(T_int2);
[dd_x3, dd_t3] = dataPros.getTimeIntervalData(T_int3);
generateStateDataPlot("State Data",dd_x1,dd_t1)
generateStateDataPlot("State Data",dd_x2,dd_t2)
generateStateDataPlot("State Data",dd_x3,dd_t3)

%% Save Data-driven data
load("data/X_T_2407171817.mat")
X{1} = dd_x1;
X{2} = dd_x2;
X{3} = dd_x3;
% X{4} = dd_x1;
% X{5} = dd_x2;
% X{6} = dd_x3;
% X{7} = dd_x1;
% X{8} = dd_x2;
% X{9} = dd_x3;
% X{10} = dd_x1;
% X{11} = dd_x2;
% X{12} = dd_x3;
%% 
% filename = "2407171817";
% save("data/X_T_" + filename + ".mat", "X", "T", "K", "STEP_SIZE")
%%
X = cell(12,1);
T = cell(12,1);
T{1} = dd_t1;
T{2} = dd_t2;
T{3} = dd_t3;
T{4} = dd_t1;
T{5} = dd_t2;
T{6} = dd_t3;
T{7} = dd_t1;
T{8} = dd_t2;
T{9} = dd_t3;
T{10} = dd_t1;
T{11} = dd_t2;
T{12} = dd_t3;

%% Plot ALL Data
% figure("Name","All Data")
% sgtitle("All Data",'Interpreter','latex')
% for i = 1:12
%     subplot(4,3,i)
%     plot(T{i},X{i},'.-',MarkerSize=10)
%     xlabel('$t$','Interpreter','latex')
%     ylabel('$x$','Interpreter','latex')
%     xlim([T{i}(1),T{i}(end)])
% end

%% Get Data-driven State Data
% idx_switchON = find(switch_state.Data(1,:),1); % get switch on time
% % l = 10; % segment(equation) numbers
% % T_int = 0.05;   % integral time
% % N_data_points = l * T_int / STEP_SIZE;
% T_data = 1.2;
% N_data_points = T_data / STEP_SIZE;
% 
% dd_time = time.Data(idx_switchON:idx_switchON + N_data_points)';
% dd_switch_state = switch_state.Data(1,idx_switchON:idx_switchON + N_data_points);
% dd_x = state4dim.Data(idx_switchON:idx_switchON + N_data_points,:)';
% % for old data (before 240711)
% % dd_x(1,:) = motor_angle.Data(time_on_idx:time_on_idx + N_data_points);
% % dd_x(2,:) = motor_speed.Data(time_on_idx:time_on_idx + N_data_points);
% % dd_x(3,:) = pendulum_angle.Data(time_on_idx:time_on_idx + N_data_points);
% % dd_x(4,:) = pendulum_speed.Data(time_on_idx:time_on_idx + N_data_points);
% 
% 
% % show data to compute data equation
% generatePhyiscalDataPlot("Data-driven Data",dd_x,dd_time,dd_switch_state)
% 
% generateStateDataPlot("State Data",dd_x,dd_time)

%% Data Smoothing 
% dd_x = smoothdata(dd_x,2,"gaussian",20);
% smooth_mehod = "gaussian";
% winsize = 10;
% dd_x(2,:) = smoothdata(dd_x(2,:),num2str(smooth_mehod),winsize);
% dd_x(4,:) = smoothdata(dd_x(4,:),num2str(smooth_mehod),winsize);

%% Generate Plots Fcns
function generatePhyiscalDataPlot(title,x,t,swt)
figure("Name",title)
sgtitle(title,'Interpreter','latex')
subplot(3,1,1)
hold on
plot(t,x(1,:),'.-',MarkerSize=10)
plot(t,x(3,:),'.-',MarkerSize=10)
legend("motor angle","pendulum angle", ...
    'Interpreter','latex',Location="northwest")
xlim([t(1),t(end)])
subplot(3,1,2)
hold on
plot(t,x(2,:),'.-',MarkerSize=10)
plot(t,x(4,:),'.-',MarkerSize=10)
legend("motor speed","pendulum speed",...
    'Interpreter','latex',Location="northwest")
xlim([t(1),t(end)])
subplot(3,1,3)
hold on
plot(t,swt,'.-',MarkerSize=10)
legend("switch state",...
    'Interpreter','latex',Location="northwest")
xlim([t(1),t(end)])
end

function generateStateDataPlot(title,x,t)
figure("Name",title)
sgtitle(title,'Interpreter','latex')
plot(t,x,'.-',MarkerSize=10);
xlabel('$t$','Interpreter','latex')
ylabel('$x$','Interpreter','latex')
xlim([t(1),t(end)])
ax = gca;
ax.FontSize = 14;
end