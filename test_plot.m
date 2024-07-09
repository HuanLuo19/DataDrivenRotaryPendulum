x = 1:100;
rng(0,"twister")
A = cos(2*pi*0.05*x+2*pi*rand) + 0.5*randn(1,100);

[B,winsize] = smoothdata(A,"gaussian");
winsize

C = smoothdata(A,"gaussian",20);

figure
hold on
plot(x,A)
plot(x,B)
plot(x,C)
legend("Raw Data","Small Window","Large Window")