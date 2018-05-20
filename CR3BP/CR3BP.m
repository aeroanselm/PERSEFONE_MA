clear
close all
clc

addpath(genpath('./functions'));
options = odeset('Reltol',1e-13,'AbsTol',1e-14);
mu = 0.0121506683;
y0 = [1.062 0 0 0 0.461 0 1.86];
T = 1.86;

%[t,state]=ode113(@(t,x)dyn_CR3BP1(t,x,mu),[0 T],x0,options);

x = fsolve(@(y)correction(mu,y),y0,options);

[t,state]=ode113(@(t,x)dyn_CR3BP1(t,x,mu),[0 5*x(7)],x,options);


figure()
hold on
grid on
axis equal
comet3(state(:,1),state(:,2),state(:,3))