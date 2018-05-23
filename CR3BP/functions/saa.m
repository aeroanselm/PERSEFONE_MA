function [theta,S,T] = saa(mu,r0,v0,tof)

s0 = [];
s0(1) = r0(1);
s0(2) = r0(2);
s0(3) = r0(3);
s0(4) = v0(1);
s0(5) = v0(2);
s0(6) = v0(3);

options=odeset('Reltol',1e-13,'AbsTol',1e-14);
[T,S]=ode113(@(t,y) dyn_2BP(t,y,mu),[0,tof],s0,options);

[n,m] = size(S);

theta = [];
for k = 1:n
    R = [S(k,1) S(k,2) S(k,3)];
    V = [S(k,4) S(k,5) S(k,6)];
    theta(k) = 90-acos(dot(R,V)/(norm(R)*norm(V)))*180/pi;
end

% figure()
% hold on 
% grid on
% axis equal
% plot3(S(:,1), S(:,2), S(:,3))