clear 
close all
clc

s0 = [10 -2.5 0 0.01 0.01 0];
R = 6378+400;
mu = 398600;
tau = 13*3600;
DV = [];
for i = 1:tau
    [dv] = evalrelmot(s0,i,R,mu);
    DV = [DV dv];
end

figure ()
hold on
grid on
plot(1:tau, DV)
[DV,Dv1,Dv2,s] = relmot(s0,R,12495,mu)

figure()
hold on
grid on
comet3(s(:,1),s(:,2),s(:,3))
