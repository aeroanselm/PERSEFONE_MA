function [S_d] = dynCR3BPpert_err(t,S,mu,flag)

x = S(1);
y = S(2);
z = S(3);
Vx = S(4);
Vy = S(5);
Vz = S(6);

DU = 9379.11746340512;
TU = 4389.11709709003;

r1 = sqrt((x + mu)^2 + y^2 + z^2);
r2 = sqrt((x - (1 - mu))^2 + y^2 + z^2);

[n,m] = size(S);
S_d = zeros(n,m);

J2p = 0.105;
J2m = 0.00196;

J3p = 0.00775;
J3m = 0.0000315;

J4p = 0.0229;
J4m = 0.0000154;


Rp = astroConstants(32);
Rm = astroConstants(24);

Rp_a = Rp/DU;
Rm_a = Rm/DU;

mup = astroConstants(33);
mum = astroConstants(14);

mup_a = mup/DU^3*TU^2;
mum_a = mum/DU^3*TU^2;

rp = sqrt((x-(1-mu))^2 + y^2 + z^2);
rm = sqrt((x + mu)^2 + y^2 +z^2);

% c2p = J2p*mup_a*Rp_a^2/2; 
% c2m = J2m*mum_a*Rm_a^2/2;
% 
% acc_J2p_x =  -3*c2p*(x-(1-mu))/rp^5*(1 - 5*z^2/rp^2);
% acc_J2p_y =  -3*c2p*y/rp^5*(1 - 5*z^2/rp^2);
% acc_J2p_z =  -3*c2p*z/rp^5*(3 - 5*z^2/rp^2);
% 
% acc_J2m_x =  -3*c2m*(x-(-mu))/rm^5*(1 - 5*z^2/rm^2);
% acc_J2m_y =  -3*c2m*y/rm^5*(1 - 5*z^2/rm^2);
% acc_J2m_z =  -3*c2m*z/rm^5*(3 - 5*z^2/rm^2);

c2p = -3/2*J2p*(mup_a/rp^2)*(Rp_a/rp)^2;
c2m = -3/2*J2m*(mum_a/rm^2)*(Rm_a/rm)^2;

acc_J2p_x = c2p*((1 - 5*(z/rp)^2)*(x-(1-mu))/rp);
acc_J2p_y = c2p*((1 - 5*(z/rp)^2)*y/rp);
acc_J2p_z = c2p*((3 - 5*(z/rp)^2)*z/rp);

acc_J2m_x = c2m*((1 - 5*(z/rm)^2)*(x + mu)/rm);
acc_J2m_y = c2m*((1 - 5*(z/rm)^2)*y/rm);
acc_J2m_z = c2m*((3 - 5*(z/rm)^2)*z/rm);

c3p = 1/2*J3p*(mup_a/rp^2)*(Rp_a/rp)^3;
c3m = 1/2*J3m*(mum_a/rm^2)*(Rm_a/rm)^3;

acc_J3p_x = c3p*(5*(7*(z/rp)^3 - 3*(z/rp))*(x-(1-mu))/rp);
acc_J3p_y = c3p*(5*(7*(z/rp)^3 - 3*(z/rp))*y/rp);
acc_J3p_z = c3p*(3*(1 - 10*(z/rp)^2 + 35/3*(z/rp)^4));

acc_J3m_x = c3m*(5*(7*(z/rm)^3 - 3*(z/rm))*(x + mu)/rm);
acc_J3m_y = c3m*(5*(7*(z/rm)^3 - 3*(z/rm))*y/rm);
acc_J3m_z = c3m*(3*(1 - 10*(z/rm)^2 + 35/3*(z/rm)^4));

c4p = 5/8*J4p*(mup_a/rp^2)*(Rp_a/rp)^4;
c4m = 5/8*J4m*(mum_a/rm^2)*(Rm_a/rm)^4;

acc_J4p_x = c4p*((3 - 42*(z/rp)^2 + 63*(z/rp)^4)*(x-(1-mu))/rp);
acc_J4p_y = c4p*((3 - 42*(z/rp)^2 + 63*(z/rp)^4)*y/rp);
acc_J4p_z = c4p*((15 - 70*(z/rp)^2 + 63*(z/rp)^4)*z/rp);

acc_J4m_x = c4m*((3 - 42*(z/rm)^2 + 63*(z/rm)^4)*(x + mu)/rm);
acc_J4m_y = c4m*((3 - 42*(z/rm)^2 + 63*(z/rm)^4)*y/rm);
acc_J4m_z = c4m*((15 - 70*(z/rm)^2 + 63*(z/rm)^4)*z/rm);

S_d(1) = S(4);
S_d(2) = S(5);
S_d(3) = S(6);
S_d(4) = x + 2*Vy - (1 - mu)/r1^3*(x + mu) - mu/r2^3*(x - (1 - mu))...
    + acc_J2p_x + acc_J2m_x + acc_J3p_x + acc_J3m_x + acc_J4p_x + acc_J4m_x;
S_d(5) = y -2*Vx - y*((1 - mu)/r1^3 + mu/r2^3)...
    + acc_J2p_y + acc_J2m_y + acc_J3p_y + acc_J3m_y + acc_J4p_y + acc_J4m_y;
S_d(6) = -z*((1 - mu)/r1^3 + mu/r2^3)...
    + acc_J2p_z + acc_J2m_z + acc_J3p_z + acc_J3m_z + acc_J4p_z + acc_J4m_z;

if flag == 1
    plot(t,acc_J2p_x,'or');
    plot(t,acc_J3p_x,'og');
    plot(t,acc_J4p_x,'ob');
end
return