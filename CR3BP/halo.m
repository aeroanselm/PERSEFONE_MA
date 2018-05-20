clear
close all
clc

options = odeset('Reltol',1e-13,'AbsTol',1e-14);
mu = 0.012150586550569;
y0=[0.862 0 0.185 0 0.25 0 1.132];

x = fsolve(@(y)correctionhalo(mu,y),y0);
[t,state]=ode113(@(t,x)dyn_CR3BP1(t,x,mu),[0 5*x(7)],x,options);


figure()
hold on
grid on
axis equal
plot3(state(:,1),state(:,2),state(:,3))
