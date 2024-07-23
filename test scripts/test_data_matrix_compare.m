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


