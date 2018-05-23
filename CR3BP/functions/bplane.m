function [R,T,S] = bplane(vv_inf)

S = vv_inf/norm(vv_inf);
k = [0, 0, 1];
T = cross(S,k)/norm(cross(S,k));
R = cross(S,T);
return                 