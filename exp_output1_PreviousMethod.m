clear variables
close all
clc
%%
date = "250202";
filename = "RawData_250202";
step = 1;
figure
hold on

sum_K2_repeat = [];
for repeat = 1:10
    % repeat = 1; % change segment_number below if [step,repeat] = [4,2] [7,4] [7,5] [8,4] [8,5]
    
    load("data\data_" + date + "\Step_" + num2str(step) + "\" + ...
        filename + "_" + num2str(step) + "." + num2str (repeat) + ".mat")
    
    %% Get Data-driven Dataset
    % plot raw data: state, switch
    
    legend_name_experiment = ["motor angle","motor speed","pendulum angle","pendulum speed"];
    
    % divide segments
    idx_start = [];
    idx_end = [];
    k_start = 1;
    k_end = 1;
    for i = 1:size(time.data,1) - 1 
        if switch_state.Data(1,i) == 0 && switch_state.Data(1,i+1) == 1
            idx_start(k_start) = i + 1;
            k_start = k_start + 1;
        end
        if switch_state.Data(1,i) == 1 && switch_state.Data(1,i+1) == 0
            idx_end(k_end) = i;
            k_end = k_end + 1;
        end
    end
    segment_number = size(idx_start,2);
    % segment_number = 20; % for step 4, repeat 2
    % segment_number = 19; % for step 7, repeat 4
    % segment_number = 13; % for step 7, repeat 5
    % segment_number = 14; % for step 8, repeat 4
    % segment_number = 19; % for step 8, repeat 5
    
    % construct raw data set
    X_raw = cell(segment_number,1);
    T_raw = cell(segment_number,1);
    for i = 1:segment_number
        x_raw = state4dim.Data(idx_start(i):idx_end(i),:)';
        t_raw = time.Data(idx_start(i):idx_end(i)) - time.Data(idx_start(i));
        X_raw{i} = x_raw;
        T_raw{i} = t_raw;
    end
    % figure("Name","State segments - Raw")
    % sgtitle("State segments - Raw",'Interpreter','latex')
    % for i = 1:segment_number
    %     subplot(4,5,i)
    %     plot(T_raw{i},X_raw{i},'.-',MarkerSize=10)
    %     xlabel('$t$','Interpreter','latex')
    %     ylabel('$x$','Interpreter','latex')
    %     xlim([T_raw{i}(1),T_raw{i}(end)])
    %     grid on
    % end
    
    % get dd data: set time interval from raw data set
    T_int = [0.5,2.5];
    X = cell(segment_number,1);
    T = cell(segment_number,1);
    for i = 1:segment_number
        dd_x = X_raw{i}(:,T_int(1)/STEP_SIZE:T_int(2)/STEP_SIZE);
        dd_t = T_raw{i}(T_int(1)/STEP_SIZE:T_int(2)/STEP_SIZE);
        X{i} = dd_x;
        T{i} = dd_t;
    end
    % figure("Name","State segments")
    % % sgtitle("State segments",'Interpreter','latex')
    % for i = 1:segment_number
    %     subplot(4,5,i)
    %     plot(T{i},X{i},'.-',MarkerSize=8)
    %     xlabel('$t$','Interpreter','latex')
    %     ylabel('$x$','Interpreter','latex')
    %     xlim([T{i}(1),T{i}(end)])
    %     grid on
    % end
    
    %% Data-driven Solution
    % system
    A = [0    1.0000         0         0;
         0  -12.2135         0         0;
         0         0         0    1.0000;
         0   -7.6602  -66.8782   -0.2289];
    % B = [0   39.2743         0   24.6327]';
    % Q = diag([100 0 1000 0]);
    % R = 100;
    
    if ~all(eig(A - B * K)<0)
        fprintf("K is NOT a stabilizing gain! \n")
        return
    end 
    
    x0 = zeros(size(A,1),1);
    sys = linearSys(A,B,x0,Q,R);
    
    % ddlyap solution
    l = 0;
    epsilon = 1;
    eps_dd = eps_ddLyap(X, l, STEP_SIZE, P0, K, sys, epsilon);
    
    % ---Pi---
    eps_Pi = eps_dd.Pi_e;
    eps_Pi_lyap = eps_dd.Pi_lyap_e;
    eps_Delta_Pi = eps_Pi - eps_Pi_lyap;
    eps_norm_Delta_Pi = norm(eps_Delta_Pi);
    % fprintf(['|eps_Delta_Pi| = ', '\n'])
    % fprintf([num2str(eps_norm_Delta_Pi), '\n'])
    
    % ---Ki---
    K
    eps_Kip1 = eps_dd.Kip1_e
    eps_Kip1_lyap = eps_dd.Kip1_lyap_e
    K_opt = sys.K_opt;
    eps_Delta_K = K - K_opt;
    eps_norm_Delta_K = norm(eps_Delta_K);
    fprintf(['|K_i - K^*| = ', '\n'])
    fprintf([num2str(eps_norm_Delta_K), '\n'])
    
    eps_Delta_Ki = eps_Kip1 - K_opt;
    eps_norm_Delta_Ki = norm(eps_Delta_Ki);
    fprintf(['|K_i+1,r - K^*| = ', '\n'])
    fprintf([num2str(eps_norm_Delta_Ki), '\n'])
    
    %%
    
    sum_K2_repeat(repeat,:) = eps_Kip1;
    % plot([1 step+1],[norm(K - K_opt) norm(eps_Kip1 - K_opt)],'bx-',MarkerSize=12,LineWidth=1.2)
end
K2_average = sum(sum_K2_repeat,1)/10;
%%
figure
hold on
for repeat = 1:10
    % h2 = plot([1 step+1],[norm(K - K_opt) norm(sum_K2_repeat(repeat,:) - K_opt)],'bx',MarkerSize=12,LineWidth=1.2);
    h1 = plot(step+1,norm(sum_K2_repeat(repeat,:) - K_opt),'bx',MarkerSize=12,LineWidth=1.2);
end
h2 = plot([1 step+1], [norm(K - K_opt) norm(K2_average - K_opt)],'g.-',MarkerSize=26,LineWidth=1.2);
xlabel("$i$","Interpreter","latex")
legend([h1,h2],"$\| K_i - K^* \| $ Candidates","$\| K_i - K^* \|$ Average", 'Location','northeast', "Interpreter","latex")
xlim([1 9])
ylim([0 1000])
grid on
box on
ax = gca;
% ax.FontSize = 14;
ax.XTick = 1: 10;