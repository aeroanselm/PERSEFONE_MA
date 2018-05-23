clear
close all
clc
A = []; b = []; Aeq = []; beq = [];
lb = 1*3600;
ub = 5*3600;
tau0 = 5*3600;
options = optimoptions('fmincon','ConstraintTolerance',1e-3,...
                        'PlotFcn',{@optimplotx,@optimplotfunccount,@optimplotfval});
tau = fmincon(@(tau)relmot([0.1 0.02 0 0 0 0],6378+400,tau0,398600),tau0,...
                      A,b,Aeq,beq,lb,ub,options);

dv=[];                  
for i = lb:ub                  
[DV,Dv1,Dv2] = relmot([0.1 0.44 0 -2e-3 -6e-3 0],6378+400,i,398600);
dv = [dv DV];
end

figure()
hold on
grid on
plot(lb:ub,dv)

[DV,Dv1,Dv2,s] = relmot([4 1 0 2e-3 6e-3 0],6378+400,3600,398600);
figure()
hold on
grid on
axis equal
plot3(s(:,1),s(:,2),s(:,3))

