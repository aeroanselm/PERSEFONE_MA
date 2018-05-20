clear
close all
clc

load('./data/arrival.mat');
load('./data/dateesc.mat');
load('./data/sLMO.mat');
load('./data/circ.mat');
dt = DAYa:date_escape;

[~,m] = size(dt);

beta = zeros(1,m);
for i = 1:m
    beta(1,i) = betaangle(s5(1,:),dt(i));
end
    
figure()
hold on
grid on
plot(dt,beta*180/pi,'Linewidth',2);
xlabel('Date [mjd2000]','Interpreter','latex');
ylabel('$\beta$ angle [deg]','Interpreter','latex');
title('\textbf{$\beta$ angle vs date: Phobos orbit}','Interpreter','latex');

dt = date_escape:(date_escape + 4*365);
[~,m] = size(dt);

beta = zeros(1,m);
for i = 1:m
    beta(1,i) = betaangle(statelmo(1,:),dt(i));
end
    
figure()
hold on
grid on
plot(dt,beta*180/pi,'Linewidth',2);
xlabel('Date [mjd2000]','Interpreter','latex');
ylabel('$\beta$ angle [deg]','Interpreter','latex');
title('\textbf{$\beta$ angle vs date: LMO}','Interpreter','latex');
