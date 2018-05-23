function [c, ceq] = mancon(x,T,R,rsoi,vvsoi,mu,kep_p,date_soi)

Delta = x(1);
theta = x(2);

DeltaT = Delta*cos(theta)*T;
DeltaR = Delta*sin(theta)*R;
Delta = DeltaT+DeltaR;
[DVs, DVtot, T, Ttot, rp_H] = manoeuvres(Delta, rsoi, vvsoi, mu, kep_p, date_soi);

ceq = 0;
c = (astroConstants(24)+120)-rp_H;