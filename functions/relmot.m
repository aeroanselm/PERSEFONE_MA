function [DV,Dv1,Dv2,s] = relmot(s0,rt,tau,mu)

x0   = [s0(1) s0(2) s0(3)];
x0_d = [s0(4) s0(5) s0(6)];
om = sqrt(mu/rt^3);

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
DV = Dv1+Dv2;

s = [];
t = 0:1:tau;
for k = 1:length(t)
    s(k,1) = x0(1)+2*x0r_d(2)/om*(1-cos(om*t(k)))+(4*x0r_d(1)/om-6*x0(2))*sin(om*t(k))+(6*om*x0(2)-3*x0r_d(1))*t(k);
    s(k,2) = 4*x0(2)-2*x0r_d(1)/om+(2*x0r_d(1)/om-3*x0(2))*cos(om*t(k))+x0r_d(2)/om*sin(om*t(k));
    s(k,3) = x0r_d(3)*cos(om*t(k))-x0(3)*om*sin(om*t(k));
end






