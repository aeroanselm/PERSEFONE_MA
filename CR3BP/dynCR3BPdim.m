function [sd] = dynCR3BPdim (t,s,d,m1,m2,mu1,mu2)

G = astroConstants(1);                  
M = m1 + m2;
mu = G*M;
OM = sqrt(mu/d^3);
p1 = m1/M;
p2 = m2/M;

x1 = -p2*d;
x2 = p1*d;
x = s(1);
y = s(2);
z = s(3);
Vx = s(4);
Vy = s(5);
Vz = s(6);

r1 = sqrt((x-x1)^2 + y^2 + z^2);
r2 = sqrt((x-x2)^2 + y^2 + z^2);
[n,m] = size(s);
sd = zeros(n,m);
sd(1) = Vx;
sd(2) = Vy;
sd(3) = Vz;
sd(4) = -mu1/r1^3*(x + p2*d) - mu2/r2^3*(x - p1*d) + 2*OM*Vy + OM^2*x;
sd(5) = -mu1/r1^3*y - mu2/r2^3*y - 2*OM*Vx + OM^2*y;
sd(6) = -mu1/r1^3*z - mu2/r2^3*z;
return