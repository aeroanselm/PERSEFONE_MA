function [S_d] = dynCR3BPpert(t,S,mu)

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

Rp = astroConstants(32);
Rm = astroConstants(24);

Rp_a = Rp/DU;
Rm_a = Rm/DU;

mup = astroConstants(33);
mum = astroConstants(14);

mup_a = mup/DU^3*TU^2;
mum_a = mum/DU^3*TU^2;

c1 = J2p*mup_a*Rp_a^2/2; 
c2 = J2m*mum_a*Rm_a^2/2;

rp = sqrt((x-(1-mu))^2 + y^2 + z^2);
rm = sqrt((x + mu)^2 + y^2 +z^2);

acc_J2p_x =  3*c1*(x-(1-mu))/rp^5*(1 - 5*z^2/rp^2);
acc_J2p_y =  3*c1*y/rp^5*(1 - 5*z^2/rp^2);
acc_J2p_z =  3*c1*z/rp^5*(3 - 5*z^2/rp^2);

acc_J2m_x =  3*c2*(x-(-mu))/rm^5*(1 - 5*z^2/rm^2);
acc_J2m_y =  3*c2*y/rm^5*(1 - 5*z^2/rm^2);
acc_J2m_z =  3*c2*z/rm^5*(3 - 5*z^2/rm^2);


S_d(1) = S(4);
S_d(2) = S(5);
S_d(3) = S(6);
S_d(4) = x + 2*Vy - (1 - mu)/r1^3*(x + mu) - mu/r2^3*(x - (1 - mu)) + acc_J2p_x + acc_J2m_x;
S_d(5) = y -2*Vx - y*((1 - mu)/r1^3 + mu/r2^3) + acc_J2p_y + acc_J2m_y;
S_d(6) = -z*((1 - mu)/r1^3 + mu/r2^3) + acc_J2p_z + acc_J2m_z;


return