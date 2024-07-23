%% Load Data 
clear variables
close all
clc
% --- load raw data
load("data\data_240716\data_raw_2407161541.mat")
% --- Show Raw Data Plot
generatePhyiscalDataPlot("Raw Data",state4dim.Data',time.Data,switch_state.Data(1,:))

%% Choose Animation Data
chooseAnimData = dataProcessing(state4dim.Data',time.Data',switch_state.Data(1,:));
% Choose animation time
T_int = [5,20]; % [start, end] time from raw data time
[anim_x_raw, anim_t_raw, anim_swt_raw] = chooseAnimData.getTimeIntervalData(T_int);
% change data to 1/25 step size (frame rate)
framerate = 25;
idx_25fr = (T_int(2) - T_int(1))/STEP_SIZE / (1/framerate/STEP_SIZE);
anim_x_25fr = zeros(4,idx_25fr);
anim_t_25fr = zeros(4,idx_25fr);
anim_swt_25fr = false(1,idx_25fr);
for i = 1:idx_25fr
    anim_x_25fr(:,i) = anim_x_raw(:,1 + i*(1/framerate/STEP_SIZE)); % data setp size = 0.005s, frame rate = 1/25s
    anim_t_25fr(:,i) = anim_t_raw(:,1 + i*(1/framerate/STEP_SIZE));
    anim_swt_25fr(i) = anim_swt_raw(1 + i*(1/framerate/STEP_SIZE));
end
%% Generate Animation
generatePendulumAnim(anim_x_25fr(1,:)',anim_x_25fr(3,:)',anim_swt_25fr,anim_t_25fr,"animation\test")

%% FUNCTIONS
% Generate Animation Fcn
function generatePendulumAnim(motor_angle, pendulum_angle, switch_state, time, filename)
videoObject = VideoWriter(filename,'MPEG-4');
videoObject.FrameRate = 25;
open(videoObject)

base_height = 3;
arm_length = 1;
pendulem_length = 2;
a = base_height;
b = arm_length;
c = pendulem_length;

th1 = motor_angle;
phi2 = pendulum_angle;

A = [0, 0, a];
B = [b .* cos(th1), ...
    b .* sin(th1), ...
    a .* ones(size(motor_angle))];
C = [b .* cos(th1) - c .* sin(phi2) .* sin(th1), ...
    b .* sin(th1) - c .* sin(phi2) .* cos(th1), ...
    a.* ones(size(motor_angle)) - c .* cos(phi2)];

figure
[X,Y,Z] = sphere;
handle_sphere = surf(X/2-1,Y/2+1,Z/2+2);
shading interp
handle_sphere.FaceColor = "g";
hold on
plot3([0; A(1)], [0; A(2)], [0; A(3)],'.-',LineWidth=3);
handle_l1 = plot3([A(1); B(1,1)], [A(2); B(1,2)], [A(3); B(1,3)],'.-', ...
    LineWidth=3);
handle_l2 = plot3([B(1,1); C(1,1)], [B(1,2); C(1,2)], [B(1,3); C(1,3)],'.-', ...
    LineWidth=3);
xlim([-2 2])
ylim([-2 2])
zlim([0 4])
grid on
daspect([1 1 1])

for i = 1:length(motor_angle)
    set(handle_l1,'XData',[0; B(i,1)]);
    set(handle_l1,'YData',[0; B(i,2)]);
    set(handle_l1,'ZData',[A(3); B(i,3)]);

    set(handle_l2,'XData',[B(i,1); C(i,1)]);
    set(handle_l2,'YData',[B(i,2); C(i,2)]);
    set(handle_l2,'ZData',[B(i,3); C(i,3)]);
    if switch_state(i) == 1
        set(handle_sphere,'FaceColor','r')
    else
        set(handle_sphere,'FaceColor','g')
    end
    drawnow;
    writeVideo(videoObject, getframe);
end
close(videoObject)
end

% Generate Plots Fcns
function generatePhyiscalDataPlot(title,x,t,swt)
figure("Name",title)
sgtitle(title,'Interpreter','latex')
subplot(3,1,1)
hold on
plot(t,x(1,:),'.-',MarkerSize=10)
plot(t,x(3,:),'.-',MarkerSize=10)
legend("motor angle","pendulum angle", ...
    'Interpreter','latex',Location="best")
xlim([t(1),t(end)])
grid on

subplot(3,1,2)
hold on
plot(t,x(2,:),'.-',MarkerSize=10)
plot(t,x(4,:),'.-',MarkerSize=10)
legend("motor speed","pendulum speed",...
    'Interpreter','latex',Location="best")
xlim([t(1),t(end)])
grid on

subplot(3,1,3)
hold on
plot(t,swt,'.-',MarkerSize=10)
xlabel('$t$','Interpreter','latex')
legend("switch state",...
    'Interpreter','latex',Location="best")
xlim([t(1),t(end)])
end