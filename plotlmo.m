clear all 
close all
clc
%%
options = odeset('Reltol',1e-13,'AbsTol',1e-14);

mu = astroConstants(14);           % Mars' planetary constant
T_phobos = period(kep2esc,mu);
statecp = kep2car(kep2esc,mu);
[~,scp]=ode113(@(t,y)dyn_2BP(t,y,mu),[0 T_phobos],statecp,options);
Tellisse = period(kep_ellisse,mu);
state_ellisse = kep2car(kep_ellisse,mu);
[~,ellisse]=ode113(@(t,y)dyn_2BP(t,y,mu),[0 Tellisse],state_ellisse,options);
statelmo = kep2car(keplmo,mu);
Tlmo = period(keplmo,mu);
[~,lmo]=ode113(@(t,y)dyn_2BP(t,y,mu),[0 Tlmo],statelmo,options);
kep_ellisse(6)=0;
s1 = kep2car(kep_ellisse,mu);
vv1 = s1(4:6);
vu = vv1/norm(vv1);
vv2 = 3.7*vu;
s2 = [s1(1:3) vv2];
[~,ellisse2]=ode113(@(t,y)dyn_2BP(t,y,mu),[0 Tellisse],s2,options);
vv3 = 3.5*vu;
s3 = [s1(1:3) vv3];
[~,ellisse3]=ode113(@(t,y)dyn_2BP(t,y,mu),[0 Tellisse],s3,options);

figure()
hold on
axis equal
title('\textbf{LMO insertion}','Interpreter','latex');
drawPlanet('Mars',[0 0 0],gca,1);
q1 = plot3(scp(:,1),scp(:,2), scp(:,3),'Linewidth',1.5);
q2 = plot3(ellisse(:,1),ellisse(:,2), ellisse(:,3),'Linewidth',1.5);
q3 = plot3(ellisse2(:,1),ellisse2(:,2), ellisse2(:,3),'Linewidth',1.5);
q4 = plot3(ellisse3(:,1),ellisse3(:,2), ellisse3(:,3),'Linewidth',1.5);
q5 = plot3(lmo(:,1),lmo(:,2), lmo(:,3),'Linewidth',1.5);
legend([q1 q2 q3 q4 q5],'Location','southeast','Starting orbit','First transfer ellipse','First break','Second break','Final Orbit [LMO]');

