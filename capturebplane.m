function [Delta] = capturebplane(VVsoi,RRsoi,rp,mu)

rsoi = norm(RRsoi);
vsoi = norm(VVsoi);

kep = car2kep([RRsoi VVsoi],mu);
a = kep(1);
e = kep(2);
i = kep(3);
Om = kep(4);
w = kep(5);
ta = kep(6);

% ROTATIONAL MATRIX

%[Om]_k
R3 = [cos(Om), sin(Om), 0;
     -sin(Om), cos(Om), 0;
      0,      0,        1];
%[i]_i
R2 = [1     ,0        ,0;
      0     ,cos(i)   ,sin(i);
      0     ,-sin(i)  ,cos(i)];
%[w]_k
R1 = [cos(w)  ,sin(w) ,0;
      -sin(w) ,cos(w) ,0;
       0      ,0      ,1];

R = [R1*R2*R3];
Rinv = R';

% Transformation from Areocentric ecliptic to perifocal frame
Secl = [RRsoi VVsoi]';
Sper = R*Secl;



Delta = abs(kep(1))*sqrt(kep(2)^2-1);


alfa_range = [-10:0.2:10]*pi/180;






ht = cross(RRsoi,VVsoi)/norm(cross(RRsoi,VVsoi));
S = VVsoi/norm(VVsoi);
B = cross(S,ht);
k = [0 0 1];                            % Unit vector normal to the ecliptic plane
T = cross(S,k)/norm(cross(S,k));
R = cross(S,T);
alfa = asin(dot(B,R));