function [beta] = betaangle(state,date)

mu = astroConstants(4);
rr = state(1:3);
vv = state(4:6);

kep = uplanet(date,4);
SS = kep2car(kep,mu);

pp = SS(1:3)/norm(SS(1:3));
hh = cross(rr,vv)/norm(cross(rr,vv));

beta = pi/2 - acos(dot(pp,hh));
