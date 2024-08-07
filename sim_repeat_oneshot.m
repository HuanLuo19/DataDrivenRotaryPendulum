clear variables
close all
clc

%%
load('data\data_240729\X_T_2407291757.mat')

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

%%

for i = 1:1 %size(X,1)
    x0 = x0_all(:,i);
    run("sim_linearSys_oneshot.m")
end