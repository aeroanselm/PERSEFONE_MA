function [dvv, dv, dt, kep2] = changeOM(kep, mu, dOM)

a = kep(1);
e = kep(2);
i = kep(3);
OM1 = kep(4);
w1 = kep(5);
ta1 = kep(6);
if ta1 < 0
    ta1 = 2*pi +ta1;
    kep(6) = ta1;
end


alfa = acos(cos(i)^2 + sin(i)^2*cos(dOM));
u1 = asin(sin(dOM)*sin(i)/sin(alfa));
u2 = u1;
if dOM > 0
    theta1 = u1-w1;
    if theta1 < 0
        theta1 = 2*pi + theta1;
    end
    theta2 = theta1;
    w2 = u2 - theta2;
else
    theta1 =2*pi - u1 -w1;
    if theta1 < 0
        theta1 = 2*pi + theta1;
    end
    theta2 = theta1;
    w2 = 2*pi - u2 - theta2;
end

% dv calculation
P = a*(1 - e^2);
vtheta = sqrt(mu/P)*(1 + e*cos(theta1));
dv = 2*vtheta*sin(alfa/2);

% dvv calculation
OM2 = OM1+dOM;
kep1 = [a e i OM1 w1 theta1];
kep2 = [a e i OM2 w2 theta2];
state1 = kep2car(kep1,mu);
state2 = kep2car(kep2,mu);
dvv = state2(4:6)-state1(4:6);

% dt calculation
dtheta = theta1 - ta1;
if dtheta < 0
    dtheta_new = 2*pi + dtheta;
    dt = dtheta2dt(kep,mu,dtheta_new);
else
    dt = dtheta2dt(kep,mu,dtheta);
end

return


    
