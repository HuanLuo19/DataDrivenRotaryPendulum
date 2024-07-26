sim_end_point = reshape(data_end_point,[10*10,1]);
sim_int = reshape(data_int,[10*16,1]);

save("test scripts\test data matrix compare_temp",'sim_int','sim_end_point')
%%
realdata_end_point = reshape(data_end_point,[10*10,1]);
realdata_int = reshape(data_int,[10*16,1]);
load('test scripts\test data matrix compare_temp.mat')
%%

save("test scripts\test data matrix compare_2407192348",'sim_int','sim_end_point','realdata_end_point','realdata_int')
%%
load('test scripts\test data matrix compare_2407192348.mat')
figure
subplot(2,1,1)
plot(sim_end_point,'x-',MarkerSize=6)
hold on
plot(realdata_end_point,'.-',MarkerSize=6)
legend('simulation data','real data')
title('end point')
grid on

subplot(2,1,2)
plot(sim_int,'x-',MarkerSize=6)
hold on
plot(realdata_int,'.-',MarkerSize=6)
legend('simulation data','real data')
title('integral')
grid on

%% 20240725 compare linear and non-linear simulation data-driven result
n = size(X{1},1);
Qk_linear = dd_sim_linear.Q + dd_sim_linear.Ki' * dd_sim_linear.R * dd_sim_linear.Ki;
Qk_NL = dd_sim_NL.Q + dd_sim_NL.Ki' * dd_sim_NL.R * dd_sim_NL.Ki;
vec_Qk_linear = reshape(Qk_linear,[n*n,1]);
vec_Qk_NL = reshape(Qk_NL,[n*n,1]);

vecs_Pi_linear = pinv(dd_sim_linear.data_matrix_end_point) * - dd_sim_linear.data_array_integral * vec_Qk_linear;
vecs_Pi_NL = pinv(dd_sim_NL.data_matrix_end_point) * - dd_sim_NL.data_array_integral * vec_Qk_NL;
Pi_linear = vecs2mat(vecs_Pi_linear);
Pi_NL = vecs2mat(vecs_Pi_NL);

error_delta_xx = dd_sim_linear.data_matrix_end_point - dd_sim_NL.data_matrix_end_point;
error_I_xx = dd_sim_linear.data_array_integral - dd_sim_NL.data_array_integral;
error_Pi = Pi_linear - Pi_NL;
close all
figure
subplot(3,1,1)
heatmap(error_delta_xx)
title('Error of \delta_{xx,i}')
subplot(3,1,2)
heatmap(error_I_xx)
title("Error of I_{xx,i}")
subplot(3,1,3)
heatmap(error_Pi)
title("Error P_i")

%% --------------------------------------------
function [P] = vecs2mat(vecs_P)
% VECS2MAT
% P must be symmetric
% P = [a b;
%      b c]
% vecs_P = [a, 2b, c]'
n = (-1 + sqrt(1 + 8*size(vecs_P,1))) / 2;
P = zeros(n,n);
temp = 1;
for i = 1:n
    for j = i:n
        if i==j
            P(i,j) = vecs_P(temp);
        else
            P(i,j) = vecs_P(temp) / 2;
            P(j,i) = P(i,j);
        end
        temp = temp + 1;
    end
end
end
