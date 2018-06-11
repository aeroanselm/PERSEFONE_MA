function [dy] = correction_dro(mu,y0,flag)

x0 = y0(1:6);
%x0 = [y0(1) 0 0 0 y0(2) 0];
t = y0(7);
options = odeset('Reltol',1e-13,'AbsTol',1e-14);

[~,state]=ode113(@(t,x)dynCR3BPpert(t,x,mu,0),[0 t],x0,options);
x = state(end,:);

if flag == 0
    dy = norm([x(2) x(4) x(6)]);
elseif flag ==1
    x = x(2:6);
    x0 = x0(2:6);
    dy = norm([x - x0]);
end
    
