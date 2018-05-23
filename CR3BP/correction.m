function [dy] = correction(mu,y0)

x0 = y0(1:6);
t = y0(7);
options = odeset('Reltol',1e-13,'AbsTol',1e-14);

[~,state]=ode113(@(t,x)dyn_CR3BP1(t,x,mu),[0 t],x0,options);
x = state(end,:);

dy = [x(2) x(4) x(6)];


