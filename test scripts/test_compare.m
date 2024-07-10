
figure
plot(t_itr,x_itr,'r-',LineWidth=3,MarkerSize=10);
hold on

plot(t_itr,dd_x,'b-',LineWidth=3,MarkerSize=6);
xlabel('$t$','Interpreter','latex')
title('$x$','Interpreter','latex')
ax = gca;
ax.FontSize = 14;

