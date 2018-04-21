function [c,ceq] = dv1con(x)

% date1 = x(1);
% date2 = x(2);
% 
% p1 = 3;
% p2 = 4;
% mu = astroConstants(4);
% 
% kep1 = uplanet(p1);
% kep2 = uplanet(p2);
% 
% s1 = kep2car(kep1,mu);
% s2 = kep2car(kep2,mu);
[dv_d,dv_a] = evaldv (x,3,4)

ceq = 0;
c = dv_d-3.45;

