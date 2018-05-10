function [dvv, dv, dt, kep_heo, rp_h] = hyp2heo(state_h,mu,e_heo)

kep_h = car2kep(state_h,mu);

ah = kep_h(1);
eh = kep_h(2);
ih = kep_h(3);
OMh = kep_h(4);
wh = kep_h(5);
tah = kep_h(6);

kep1 = [ah eh ih OMh wh 0];
state_h = kep2car(kep1,mu);
rrp_h = state_h(1:3);
rp_h = norm(rrp_h);
vvp_h = state_h(4:6);
vvp_h_u = vvp_h/norm(vvp_h);

% Elliptical orbit
rrp_heo = rrp_h;
rp_heo = rp_h;
ra_heo = -rp_heo*(1+e_heo)/(e_heo-1);
a_heo = (rp_heo + ra_heo)/2;
P_heo = a_heo*(1 - e_heo^2);
vp_heo = sqrt(mu/P_heo)*(1 + e_heo);
vvp_heo = vp_heo*vvp_h_u;
state_heo = [rrp_heo vvp_heo];
kep_heo = car2kep(state_heo,mu);

% dvv calculation
dvv = vvp_heo - vvp_h;

% dv calculation 
dv = norm(dvv);

% dt calculation
dtheta = 0 - tah;
if dtheta < 0
    dtheta = 2*pi - abs(dtheta);
end
dt = dtheta2dt(kep_h,mu,dtheta);

return
