function [DVs, DVtot, T, Ttot, rp_H] = manoeuvres(Delta, rsoi, vvsoi, mu, kep_p, date_soi)

DVs = [];
T = [];
options = odeset('Reltol',1e-13,'AbsTol',1e-14);

% First manoeuvre, hyperbola-->HEO
S = vvsoi/norm(vvsoi);
gamma = asin(norm(Delta)/rsoi);
rrpar = rsoi*cos(gamma)*(-S);
r_1 = rrpar+Delta;
state0 = [r_1 vvsoi];
kep1 = car2kep(state0,mu);
P1 = kep1(1)*(1-kep1(2)^2);
ta_soi = acos((P1/norm(r_1)-1)/kep1(2));
dt1 = dtheta2dt (kep1,mu,ta_soi);  T(1)=dt1;
[~,state_H]=ode113(@(t,y)dyn_2BP(t,y,mu),0:dt1,state0,options);
rrp_H = state_H(end,1:3);
rp_H = norm(rrp_H);
vvp_H = state_H(end,4:6);
vvp_H_u = vvp_H/norm(vvp_H); % Direction of the velocity at pericenter
vp_H = norm(vvp_H);

e1 = 0.9;
rrp_e1 = rrp_H;
rp_e1 = rp_H;
ra_e1 = -rp_e1*(1+e1)/(e1-1);
a_e1 = (rp_e1+ra_e1)/2;
P_e1 = a_e1*(1-e1^2);
vp_e1 = sqrt(mu/P_e1)*(1+e1);
dv_H_e1 = vp_e1-norm(vvp_H); 
dvv_H_e1 = abs(dv_H_e1)*(-vvp_H_u); DVs = [DVs; dvv_H_e1]; 
vvp_e1 = vvp_H + dvv_H_e1;
kep_e1 = car2kep([rrp_e1 vvp_e1], mu);
ra_e1 = kep_e1(1)*(1+kep_e1(2));
T_e1 = period(kep_e1,mu);
[~,state_e1]=ode113(@(t,y)dyn_2BP(t,y,mu),0:T_e1,[rrp_e1 vvp_e1],options);
hh_e1 = cross(rrp_e1, vvp_e1);
hh_e1 = hh_e1/norm(hh_e1);
kep_e1(6)=pi;
Sa1 = kep2car(kep_e1,mu);
rra_e1 = Sa1(1:3);
vva_e1 = Sa1(4:6);
vva_e1_u = vva_e1/norm(vva_e1);
va_e1 = norm(vva_e1);
dt2 = dtheta2dt (kep_e1,mu,pi);  T(2)=dt2;

% Change plane
rra_e2 = rra_e1;
ra_e2 = ra_e1;
rp_e2 = kep_p(1);
e2 = (ra_e2-rp_e2)/(ra_e2+rp_e2);
a2 = (ra_e2+rp_e2)/2;
P_e2 = a2*(1-e2^2);
va_e2 = sqrt(mu/P_e2)*(1-e2);
delta_i = kep_p(3)-kep_e1(3);
dv_cp = sqrt(va_e1^2+va_e2^2-2*va_e1*va_e2*cos(delta_i));

vvh = va_e2*sin(delta_i)*hh_e1;
vvpar = (va_e2*cos(delta_i)-va_e1)*vva_e1_u;
dvv_cp = vvh+vvpar; DVs = [DVs; dvv_cp];

% beta = pi-asin(va_e2*sin(delta_i)/dv_cp);
% alfaa = pi-beta;
% dv_par = dv_cp*cos(alfaa)*vva_e1_unit;
% dv_norm = dv_cp*sin(alfaa)*-hh_e1;
% dvv_cp= dv_par+dv_norm;
vva_e2 = vva_e1+dvv_cp;
kep_e2 = car2kep([rra_e2 vva_e2], mu);
if (kep_e2(3)-kep_e1(3))<0 && delta_i>0
    vvh = va_e2*sin(delta_i)*-hh_e1;
    vvpar = (va_e2*cos(delta_i)-va_e1)*vva_e1_u;
    dvv_cp = vvh+vvpar;
    vva_e2 = vva_e1+dvv_cp;
    kep_e2 = car2kep([rra_e2 vva_e2], mu);
end

dOM = kep_p(4)-kep_e2(4);
alfa = acos(cos(kep_e2(3))^2 + sin(kep_e2(3))^2*cos(dOM));
u1 = asin(sin(dOM)*sin(kep_e2(3))/sin(alfa));
u2 = u1;
ta1 = u1-kep_e2(5);
ta2 = ta1;
w2 = u2-ta2;
vtta = sqrt(mu/P_e2)*(1+kep_e2(2)*cos(kep_e2(6)));
dv_cOM = 2*vtta*sin(alfa/2);
kep_e2(4) = kep_p(4);


if ta1<pi 
    dthetadOM = pi+ta1;
    dt3 = dtheta2dt (kep_e2,mu,dthetadOM);  T(3)=dt3;
else
    dthetadOM = ta1-pi;
    dt3 = dtheta2dt (kep_e2,mu,dthetadOM);  T(3)=dt3;
end

% Circularize
vc = sqrt(mu/rp_e2);
vp_e2 = sqrt(mu/P_e2)*(1+e2);
dvcirc = vc-vp_e2;

dthetac = 2*pi - ta1;
dt4 = dtheta2dt (kep_e2,mu,dthetac);  T(4)=dt4;

Ttot = dt1+dt2+dt3+dt4;
DVtot = norm(dvv_H_e1)+norm(dvv_cp)+abs(dv_cOM)+abs(dvcirc);

