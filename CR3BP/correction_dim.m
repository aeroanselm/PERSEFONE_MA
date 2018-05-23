function [dy] = correction_dim(m1,m2,mu1,mu2,d,y0)

x0 = y0(1:6);
t = y0(7);

options = odeset('Reltol',1e-13,'AbsTol',1e-14);

[t,state]=ode113(@(t,x)dynCR3BPdim(t,x,d,m1,m2,mu1,mu2),[0 t],x0,options);

x = state(end,:);

dy = [x(2) x(4) x(6)]';


