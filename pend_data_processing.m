
%% Load Raw Data 
clear variables
close all
clc

%% Load All Raw Data from Folder
dir_path = "data\data_240729\";
matlist = dir(dir_path + 'data_raw' + '*.mat');
X = cell(0,0);            % comment out this line to add data
T = cell(0,0);            % comment out this line to add data
% load("data\data_temp.mat")  % comment out this line when add the first data group
for i = 1:length(matlist)
    load(dir_path + matlist(i).name)
    dataPros_Raw2SwtichON = dataProcessing(state4dim.Data',time.Data',switch_state.Data(1,:));
    [x_switchON, time_switchON, switch_state_ON] = dataPros_Raw2SwtichON.getSwitchOnData();

    %
    % legend_name_experiment = ["motor angle","motor speed","pendulum angle","pendulum speed"];
    % generatePhyiscalDataPlot("Raw Data",state4dim.Data',time.Data,switch_state.Data(1,:),legend_name_experiment)
    % generatePhyiscalDataPlot("Switch ON Data",x_switchON,time_switchON,switch_state_ON,legend_name_experiment)

    % select time interval and segment number and segment offset size
    dataPros_SwitchON2Datadriven = dataProcessing(x_switchON,time_switchON,[]);
    T_int = [0.2,0.7];
    offset_size = 0.5;
    segment_number = 2;
    for i = 1:segment_number
        [dd_x, dd_t, ~] = dataPros_SwitchON2Datadriven.getTimeIntervalData(T_int);
        X = [X; dd_x];
        T = [T; dd_t];
        T_int = T_int + offset_size;
    end
end

%% plot all segments
figure("Name","All Data")
sgtitle("All Data",'Interpreter','latex')
for i = 1:size(X,1)
    % subplot(size(X,1)/length(matlist),length(matlist),i)
    subplot(4,10,i)
    plot(T{i},X{i},'.-',MarkerSize=10)
    xlabel('$t$','Interpreter','latex')
    ylabel('$x$','Interpreter','latex')
    xlim([T{i}(1),T{i}(end)])
    grid on
end

% get segments IC
x0_all = zeros(size(X{1},1),size(X,1));
for i = 1:size(X,1)
    x0_all(:,i) = X{i}(:,i);
end
figure
plot(1:40,x0_all(:,1:40),'*-')
xlabel('Number of $x_0$','Interpreter','latex')
legend('$x_1$','$x_2$','$x_3$','$x_4$','Interpreter','latex')
%% --- Load Single Raw Data
load("data\data_240729\data_raw_2407291613.mat")
% --- Show Raw Data Plot
legend_name_experiment = ["motor angle","motor speed","pendulum angle","pendulum speed"];
generatePhyiscalDataPlot("Raw Data",state4dim.Data',time.Data,switch_state.Data(1,:),legend_name_experiment)

%% Get Switch ON Data
dataPros_Raw2SwtichON = dataProcessing(state4dim.Data',time.Data',switch_state.Data(1,:));
[x_switchON, time_switchON, switch_state_ON] = dataPros_Raw2SwtichON.getSwitchOnData();
% --- Show Switch ON Data Plot
generatePhyiscalDataPlot("Switch ON Data",x_switchON,time_switchON,switch_state_ON,legend_name_experiment)
% generateStateDataPlot("Switch ON State",x_switchON,time_switchON)
%% Select Data from Switch ON Data
X = cell(0,0);            % comment out this line to add data
T = cell(0,0);            % comment out this line to add data
% load("data\data_temp.mat")  % comment out this line when add the first data group

% select time interval and segment number and segment offset size
dataPros_SwitchON2Datadriven = dataProcessing(x_switchON,time_switchON,[]);
T_int = [0,0.05];
offset_size = 0.05;
segment_number = 10;
for i = 1:segment_number
    [dd_x, dd_t, ~] = dataPros_SwitchON2Datadriven.getTimeIntervalData(T_int);
    X = [X; dd_x];
    T = [T; dd_t];
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

%% ---> Save Data to Temperary file
% run "dd_save_data.m": Save Data to Temperary file section

%% ---> Save Data-driven Data
% until all temp data are collected
% run "dd_save_data.m": Save Data-driven Data section

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
function generatePhyiscalDataPlot(title,x,t,swt,LegendName)
    figure("Name",title)
    sgtitle(title,'Interpreter','latex')
    subplot(3,1,1)
    hold on
    plot(t,x(1,:),'.-','DisplayName',LegendName(1),MarkerSize=10)
    plot(t,x(3,:),'.-','DisplayName',LegendName(3),MarkerSize=10)
    legend('Interpreter','latex',Location="best")
    xlabel('$t$','Interpreter','latex')
    xlim([t(1),t(end)])
    grid on
    
    subplot(3,1,2)
    hold on
    plot(t,x(2,:),'.-','DisplayName',LegendName(2),MarkerSize=10)
    plot(t,x(4,:),'.-','DisplayName',LegendName(4),MarkerSize=10)
    legend('Interpreter','latex',Location="best")
    xlabel('$t$','Interpreter','latex')
    xlim([t(1),t(end)])
    grid on
    
    subplot(3,1,3)
    hold on
    plot(t,swt,'.-',MarkerSize=10)
    legend("switch state",...
        'Interpreter','latex',Location="best")
    xlabel('$t$','Interpreter','latex')
    xlim([t(1),t(end)])
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