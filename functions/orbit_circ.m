function [dvv, dv, dt, kep_c] = orbit_circ(kep, mu)

a = kep(1);
e = kep(2);
i = kep(3);
OM = kep(4);
w = kep(5);
ta1 = kep(6);

P = a*(1 - e^2);
rp = a*(1 - e);
vp = sqrt(mu/P)*(1 + e);
vc = sqrt(mu/rp);

% dv calculation
dv = vc-vp;

% dvv calculation
kep1 = [a e i OM w 0];
state1 = kep2car(kep1,mu);
vvp = state1(4:6);
vvp_u = vvp/norm(vvp);
vvc = vc*vvp_u;
dvv = vvc - vvp;
state2 = [state1(1:3) vvc];
kep_c = car2kep(state2,mu);
dv = norm(dvv);
dtheta = 2*pi - ta1;
dt = dtheta2dt(kep,mu,dtheta);
return
