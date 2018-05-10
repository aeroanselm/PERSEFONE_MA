function [dvv, dv, dt, kep2] = changeps_apo(kep1,rp2,mu,i_f)

a = kep1(1);
e = kep1(2);
i = kep1(3);
OM = kep1(4);
w = kep1(5);
ta = kep1(6);

kep_a = [a e i OM w pi];
state1 = kep2car(kep_a,mu);
di = i_f - i;

rra = state1(1:3);
vva = state1(4:6);
va = norm(vva);
vva_u = vva/norm(vva);
hh = cross(rra,vva)/norm(cross(rra,vva));

% Second ellipse calculation
ra = norm(rra);
rra2 = rra;
ra2 = ra;
a2 = (ra2 + rp2)/2;
e2 = (ra2 - rp2)/(ra2 + rp2);
P2 = a2*(1 - e2^2);
va2 = sqrt(mu/P2)*(1-e2);

% dvv calculation
v_norm = va2*sin(di)*hh;
v_par = (va2*cos(di) - va)*vva_u;
dvv = v_norm + v_par;
vva2 = vva + dvv;
state2 = [rra2 vva2];
kep2 = car2kep(state2,mu);

if (kep2(3) - kep1(3)) < 0 && di > 0
    v_norm = va2*sin(di)*-hh;
    v_par = (va2*cos(di) - va)*vva_u;   
    dvv = v_norm + v_par;
    vva2 = vva + dvv;
    state2 = [rra2 vva2];
    kep2 = car2kep(state2,mu);
    kep2(3) = i_f;
    state2 = kep2car(kep2,mu);
    vva2 = state2(4:6);
    dvv = vva2 - vva;
elseif (kep2(3) - kep1(3)) > 0 && di < 0
    v_norm = va2*sin(di)*-hh;
    v_par = (va2*cos(di) - va)*vva_u;   
    dvv = v_norm + v_par;
    vva2 = vva + dvv;
    state2 = [rra2 vva2];
    kep2 = car2kep(state2,mu);
    kep2(3) = i_f;
    state2 = kep2car(kep2,mu);
    vva2 = state2(4:6);
    dvv = vva2 - vva;
end

% dv calculation 
dv = norm(dvv);

% dt calculation
dtheta = pi - ta;
if dtheta < 0
    dtheta = 2*pi - abs(dtheta);
end
dt = dtheta2dt(kep1,mu,dtheta);

return