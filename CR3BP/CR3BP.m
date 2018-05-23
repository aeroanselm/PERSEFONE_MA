clear
close all
clc

addpath(genpath('./functions'));
options = odeset('Reltol',1e-13,'AbsTol',1e-16);
mu = 0.0121506683;
y0 = [1.062 0 0 0 0.461 0 1.86];
T = 1.86;

%[t,state]=ode113(@(t,x)dyn_CR3BP1(t,x,mu),[0 T],x0,options);
foptions = optimoptions('fsolve','FunctionTolerance',1e-13,'OptimalityTolerance',1e-14,'StepTolerance',1e-14);

x = fsolve(@(y)correction(mu,y),y0,foptions);
[t,state]=ode113(@(t,x)dyn_CR3BP1(t,x,mu),[0 3*x(7)],x,options);
dv = state(end,4:6) - state(1,4:6);
dy = state(end,1:3) - state(1,1:3);

figure()
hold on
grid on
axis equal
plot3(state(:,1),state(:,2),state(:,3))