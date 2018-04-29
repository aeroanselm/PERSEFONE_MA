vinf = 2.45;
e = 0.9;
rp = 120+3390;
mu = 42828;

ah = -mu/vinf^2;
vph = sqrt(mu*(2*ah-rp)/(rp*ah));

rae = -rp*(e+1)/(e-1);
a = (rp+rae)/2;
P = a*(1-e^2);
vpe = sqrt(mu/P)*(1+e);
vae = sqrt(mu/P)*(1-e);
dv = vpe-vph

rphobos = 9375;
rp2 = rphobos;
ra2 = rae;
e2 = (ra2-rp2)/(ra2+rp2);
a2 = (ra2+rp2)/2;
P2 = a2*(1-e2^2);
va2 = sqrt(mu/P2)*(1-e2);
vp2 = sqrt(mu/P2)*(1+e2);
vpho = sqrt(mu/rphobos);
dv_pho = vpho-vp2
dvcp = sqrt(vae^2+va2^2+2*vae*va2*cos(26*pi/180))
dvtot = abs(dvcp)+abs(dv)+abs(dv_pho)

kep = [a e 0 0 0 0];
T = period(kep,mu);
s0 = kep2car(kep,mu);
options=odeset('Reltol',1e-13,'AbsTol',1e-14);

[~,s]=ode113(@(t,y)dyn_2BP(t,y,mu),[0 T],s0,options);
figure()
hold on
grid on
axis equal
plot3(s(:,1),s(:,2),s(:,3))