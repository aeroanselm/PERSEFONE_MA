function [S_d] = dyn_CR3BP1 (t, S, mu)

x = S(1);
y = S(2);
z = S(3);
Vx = S(4);
Vy = S(5);
Vz = S(6);

r1 = sqrt((x + mu)^2 + y^2 + z^2);
r2 = sqrt((x - (1 - mu))^2 + y^2 + z^2);

[n,m] = size(S);
S_d = zeros(n,m);

S_d(1) = S(4);
S_d(2) = S(5);
S_d(3) = S(6);
S_d(4) = x + 2*Vy - (1 - mu)/r1^3*(x + mu) - mu/r2^3*(x - (1 - mu));
S_d(5) = y -2*Vx - y*((1 - mu)/r1^3 + mu/r2^3);
S_d(6) = -z*((1 - mu)/r1^3 + mu/r2^3);


return