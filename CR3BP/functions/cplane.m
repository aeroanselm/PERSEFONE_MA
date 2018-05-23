function [dv, kep2] = cplane(kep,OM2,i2,mu)

a = kep(1);
e = kep(2);
i1 = kep(3);
OM1 = kep(4);
w1 = kep(5);
%ta1 = kep(6);
di = i2 -i1;
dOM = OM2 - OM1;

if di*dOM > 0
    alfa = acos(cos(i1)*cos(i2) + sin(i1)*sin(i2)*cos(dOM));
    u1 = asin(sin(dOM)*sin(i2)/sin(alfa));
    u2 = asin(sin(dOM)*sin(i1)/sin(alfa));
    theta1 = u1 - w1;
    theta2 = theta1;
    w2 = u2 - theta2;
elseif di*dOM < 0
    alfa = acos(cos(i1)*cos(i2) + sin(i1)*sin(i2)*cos(dOM));
    u1 = asin(sin(dOM)*sin(i2)/sin(alfa));
    u2 = asin(sin(dOM)*sin(i1)/sin(alfa));
    theta1 = 2*pi - u1 - w1;
    theta2 = theta1;
    w2 = 2*pi - u2 - theta2;
end

P = a*(1 - e^2);
Vtheta = sqrt(mu/P)*(1 + e*cos(theta1));

dv = 2*Vtheta*sin(alfa/2);
kep2 = [a e i2 OM2 w2 theta2];

