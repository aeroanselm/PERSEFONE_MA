clear 
close all
clc

load('IC_DRO.mat');

x = IC(:,1);
vy = IC(:,5);

figure()
hold on;
grid on;
title('\textbf{Initial conditions of the DRO family}','Interpreter','Latex');
xlabel('\textbf{X (DU)}','Interpreter','Latex');
ylabel('\textbf{Vy (DU/TU)}','Interpreter','Latex');
plot(x,vy);

x = IC(1:10:end,1);
vy = IC(1:10:end,5);
plot(x,vy,'o');
legend('Vy vs X','10 km step');

%% 
