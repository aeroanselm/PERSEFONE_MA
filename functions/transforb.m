function [state_T, state_1, state_2, tof_sol] = transforb(dates,planet_1, planet_2, mu)

DAYd = dates(1,1);          % Selected departure date [mjd2000]
DAYa = dates(1,2);          % Selected arrival date   [mjd2000]
tof_sol = (DAYa - DAYd)*24*3600;        % Transfer time

% Setting ode options
options = odeset('Reltol',1e-13,'AbsTol',1e-14);

% Earth - State definition at departure and arrival dates
kep_1d = uplanet(DAYd,planet_1);
state_1d = kep2car(kep_1d,mu);
T_1 = period(kep_1d,mu);
[~,state_1] = ode113(@(t,y)dyn_2BP(t,y,mu),[0,T_1],state_1d,options);
kep_1a = kep_1d;
kep_1a(6) = kep_1d(6)+dt2dtheta(kep_1d,mu,tof_sol);
state_1a = kep2car(kep_1a,mu);

% Mars - State definition at departure and arrival dates
kep_2d = uplanet(DAYd,planet_2);
state_2d = kep2car(kep_2d,mu);
T_2 = period(kep_2d,mu);
[~,state_2] = ode113(@(t,y)dyn_2BP(t,y,mu),[0,T_2],state_2d,options);
kep_2a = kep_2d;
kep_2a(6) = kep_2d(6)+dt2dtheta(kep_2d,mu,tof_sol);
state_2a = kep2car(kep_2a,mu);

% Transfer arc computation
[~,~,~,~,V1,~] = lambertMR(state_1d(1:3),state_2a(1:3),tof_sol,mu,0,0,0,0);

% Transfer orbit computation
state0 = [state_1d(1:3) V1];        % Initial state vector
kep_T = car2kep(state0,mu);     % Keplerian elements of the transfer orbit
T_T = period(kep_T,mu);         % Transfer orbit period

[~,state_T] = ode113(@(t,y)dyn_2BP(t,y,mu),[0,(DAYa-DAYd)*3600*24],state0,options);

figure();
hold on
%grid on
axis equal
%set(gca,'Color','k')
drawPlanet('Sun',[0 0 0],gca,25);
drawPlanet('Mars',state_1d(1:3),gca,1500);
drawPlanet('Earth',state_2a(1:3),gca,1500);
p1=plot3(state_1(:,1),state_1(:,2),state_1(:,3),'-r','LineWidth',2);
p2=plot3(state_2(:,1),state_2(:,2),state_2(:,3),'-b','LineWidth',2);
p3=plot3(state_T(:,1),state_T(:,2),state_T(:,3),'-k','LineWidth',2);
%p4=plot3(state_Tp(:,1),state_Tp(:,2),state_Tp(:,3),'r','LineWidth',2);
p5=plot3(state_1d(1),state_1d(2),state_1d(3),'*k');
p6=plot3(state_1a(1),state_1a(2),state_1a(3),'ok');
plot3(state_2d(1),state_2d(2),state_2d(3),'*k');
plot3(state_2a(1),state_2a(2),state_2a(3),'ok');
legend([p1 p2 p3 p5 p6],'Location','southeast','Mars orbit','Earth orbit','Transfer orbit', 'Start','Finish');
title('\textbf{Mars - Earth interplanetary leg}','Interpreter','latex');
xlabel('x [Km]','Interpreter','latex');
ylabel('y [Km]','Interpreter','latex');
zlabel('z [Km]','Interpreter','latex');