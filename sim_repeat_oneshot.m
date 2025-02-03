clear variables
close all
clc

%%
load('data\data_241020\X_T_2410202156.mat')

%% plot all segments
figure("Name","All Data")
sgtitle("All Data",'Interpreter','latex')
for i = 1:size(X,1)
    % subplot(size(X,1)/length(matlist),length(matlist),i)
    subplot(2,10,i)
    plot(T{i},X{i},'.-',MarkerSize=10)
    xlabel('$t$','Interpreter','latex')
    ylabel('$x$','Interpreter','latex')
    xlim([T{i}(1),T{i}(end)])
    grid on
end

% get segments IC
x0_all = zeros(size(X{1},1),size(X,1));
for i = 1:size(X,1)
    x0_all(:,i) = X{i}(:,1);
end
figure
plot(1:20,x0_all(:,1:20),'*-')
xlabel('Number of $x_0$','Interpreter','latex')
legend('$x_1$','$x_2$','$x_3$','$x_4$','Interpreter','latex')

%%
X_sim = cell(0,0);         
T = cell(0,0);            
for i = 1:size(X,1)
    x0 = x0_all(:,i);

    % % linear oneshot
    run("sim_linearSys_oneshot.m")
    dataPros_SimLinear2Datadriven = dataProcessing(x_sim_linear,t_sim_linear,[]);
    T_int = [0,1];
    [x_dd, t_dd, ~] = dataPros_SimLinear2Datadriven.getTimeIntervalData(T_int);
    X_sim = [X_sim; x_dd];
    T = [T; t_dd];
    

    % % non-linear oneshot
    % run("sim_nonlinearSys.m")
    % dataPros_SimLinear2Datadriven = dataProcessing(sim_NL.state4dim.Data',sim_NL.tsim',[]);
    % T_int = [0,0.5];
    % [x_dd, t_dd, ~] = dataPros_SimLinear2Datadriven.getTimeIntervalData(T_int);
    % X_sim = [X_sim; x_dd];
    % T = [T; t_dd];

    % offset_size = 0.05;
    % segment_number = 1;
    % for i = 1:segment_number
    %     [x_dd, t_dd, ~] = dataPros_SimLinear2Datadriven.getTimeIntervalData(T_int);
    %     X = [X; x_dd];
    %     T = [T; t_dd];
    %     T_int = T_int + offset_size;
    % end

    % % Plot Data-driven Data
    % figure("Name","All Data-driven Data")
    % sgtitle("All Data-driven Data - Linear System",'Interpreter','latex')
    % for i = 1:size(X_sim,1)
    %     % subplot(size(X,1)/segment_number,segment_number,i)
    %     figure
    %     plot(T{i},X_sim{i},'.-',MarkerSize=10)
    %     xlabel('$t$','Interpreter','latex')
    %     ylabel('$x$','Interpreter','latex')
    %     xlim([T{i}(1),T{i}(end)])
    %     grid on
    % end
end
%%
segment_number = 10;
figure("Name","All Data-driven Data")
sgtitle("All Data-driven Data - Linear System",'Interpreter','latex')
for i = 1:size(X_sim,1)
    subplot(size(X_sim,1)/segment_number,segment_number,i)
    plot(T{i},X_sim{i},'.-',MarkerSize=10)
    xlabel('$t$','Interpreter','latex')
    ylabel('$x$','Interpreter','latex')
    xlim([T{i}(1),T{i}(end)])
    grid on
end