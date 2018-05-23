function [dv] = evalrelmot(s0,tau,R,mu)

x0   = [s0(1) s0(2) s0(3)];
x0_d = [s0(4) s0(5) s0(6)];
om = sqrt(mu/R^3);

Delta = 3*om*tau*sin(om*tau)-8*(1-cos(om*tau));

x0r_d = [];
x0r_d(1) = (x0(1)*sin(om*tau)+x0(2)*(6*om*tau*sin(om*tau)-14*(1-cos(om*tau))))/Delta*om;
x0r_d(2) = (2*x0(1)*(1-cos(om*tau))+x0(2)*(4*sin(om*tau)-3*om*tau*cos(om*tau)))/Delta*om;
x0r_d(3) = -x0(3)/tan(om*tau)*om;

Dv1 = sqrt((x0r_d(1)-x0_d(1))^2+(x0r_d(2)-x0_d(2))^2+(x0r_d(3)-x0_d(3))^2);

xtau_d(1) = 2*x0r_d(2)*sin(om*tau)+(4*x0r_d(1)-6*om*x0(2))*cos(om*tau)+6*om*x0(2)-3*x0r_d(1);
xtau_d(2) = (3*om*x0(2)-2*x0r_d(1))*sin(om*tau)+x0r_d(2)*cos(om*tau);
xtau_d(3) = x0r_d(3)*cos(om*tau)-x0(3)*om*sin(om*tau);

Dv2 = sqrt(xtau_d(1)^2+xtau_d(2)^2+xtau_d(3)^2);
dv = Dv1+Dv2;
