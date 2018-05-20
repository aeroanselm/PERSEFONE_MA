%This routine is used to produced data and plots about the Earth link
%possibility during the whole PERSEFONE mission
clear 
close all
clc

options = odeset('Reltol',1e-13,'AbsTol',1e-14);

mu_sun = astroConstants(4);
planet_1 = 3;
planet_2 = 4;
load ('./data/departure.mat');
load ('./data/arrival.mat');
dates = [DAYd, DAYa];
[state_T, state_1, state_2, tof_sol] = transforb(dates,planet_1, planet_2, mu_sun);

[tv,visibility] = evw_sun(state_1(1,:), state_T(1,:), tof_sol);

figure()
hold on;
grid on;
plot(tv/3600/24,visibility,'b','Linewidth',3);
xlabel('time of flight [days]','Interpreter','latex');
ylabel('Earth link availability [1 = yes, 0 = no]','Interpreter','latex');
title('\textbf{Earth-Mars transfer Earth link availability}','Interpreter','latex');

% return to EARTH
load('./data/medates.mat');
[state_T, state_2, state_1, tof_sol] = transforb(medates,planet_2, planet_1, mu_sun);
[tv,visibility] = evw_sun(state_1(1,:), state_T(1,:), tof_sol);

figure()
hold on;
grid on;
plot(tv/3600/24,visibility,'b','Linewidth',3);
xlabel('time of flight [days]','Interpreter','latex');
ylabel('Earth link availability [1 = yes, 0 = no]','Interpreter','latex');
title('\textbf{Mars-Earth transfer Earth link availability}','Interpreter','latex');

load('./data/iperbole.mat');
load('./data/ellisse.mat');
load('./data/cps.mat');
load('./data/cOM.mat');
load('./data/dt.mat');
load('./data/DAYc.mat');

% Man 1
kep_Ec = uplanet(date_capture,3);
s0E = kep2car(kep_Ec,mu_sun);

kep_Mc = uplanet(date_capture,4);
s0M = kep2car(kep_Mc,mu_sun);

[tv,visibility_sun,visibility_mars] = evw_sun(s0E, s1(1,:), dt(1),s0M);

figure()
hold on;
grid on;
plot(tv/3600,visibility_sun,'y','Linewidth',3);
plot(tv/3600,visibility_mars,'r','Linewidth',3);
xlabel('time of flight [hours]','Interpreter','latex');
ylabel('Earth link availability [1 = yes, 0 = no]','Interpreter','latex');
title('\textbf{Man 1 link availability}','Interpreter','latex');
legend('Location','southeast','Sun','Mars');

% Man 2
kep_Ec = uplanet((date_capture + dt(1)/3600/24) ,3);
s0E = kep2car(kep_Ec,mu_sun);

kep_Mc = uplanet((date_capture + dt(1)/3600/24),4);
s0M = kep2car(kep_Mc,mu_sun);

[tv,visibility_sun,visibility_mars] = evw_sun(s0E, s2(1,:), dt(2),s0M);

figure()
hold on;
grid on;
plot(tv/3600,visibility_sun,'y','Linewidth',3);
plot(tv/3600,visibility_mars,'r','Linewidth',3);
xlabel('time of flight [hours]','Interpreter','latex');
ylabel('Earth link availability [1 = yes, 0 = no]','Interpreter','latex');
title('\textbf{Man 2 link availability}','Interpreter','latex');
legend('Location','southeast','Sun','Mars');

% Man 3
kep_Ec = uplanet((date_capture + dt(1)/3600/24 + dt(2)/3600/24) ,3);
s0E = kep2car(kep_Ec,mu_sun);

kep_Mc = uplanet((date_capture + dt(1)/3600/24 + dt(2)/3600/24),4);
s0M = kep2car(kep_Mc,mu_sun);

[tv,visibility_sun,visibility_mars] = evw_sun(s0E, s3(1,:), dt(3),s0M);

figure()
hold on;
grid on;
plot(tv/3600,visibility_sun,'y','Linewidth',3);
plot(tv/3600,visibility_mars,'r','Linewidth',3);
xlabel('time of flight [hours]','Interpreter','latex');
ylabel('Earth link availability [1 = yes, 0 = no]','Interpreter','latex');
title('\textbf{Man 3 link availability}','Interpreter','latex');
legend('Location','southeast','Sun','Mars');

% Man 4
kep_Ec = uplanet((date_capture + dt(1)/3600/24 + dt(2)/3600/24 + dt(3)/3600/24),3);
s0E = kep2car(kep_Ec,mu_sun);

kep_Mc = uplanet((date_capture + dt(1)/3600/24 + dt(2)/3600/24 + dt(3)/3600/24),4);
s0M = kep2car(kep_Mc,mu_sun);

[tv,visibility_sun,visibility_mars] = evw_sun(s0E, s3(1,:), dt(4),s0M);

figure()
hold on;
grid on;
plot(tv/3600,visibility_sun,'y','Linewidth',3);
plot(tv/3600,visibility_mars,'r','Linewidth',3);
xlabel('time of flight [hours]','Interpreter','latex');
ylabel('Earth link availability [1 = yes, 0 = no]','Interpreter','latex');
title('\textbf{Man 4 link availability}','Interpreter','latex');
legend('Location','southeast','Sun','Mars');

%% Escape 
load('./data/dateesc');
load('./data/ellissesc.mat');
load('./data/hypesc.mat');
load('./data/sLMO.mat');
load('./data/Tellisseesc.mat');
load('./data/Tlmo.mat');
date_escape = date_escape - 2.5;

% Man esc 1
kep_Ec = uplanet(date_escape,3);
s0E = kep2car(kep_Ec,mu_sun);

kep_Mc = uplanet(date_escape,4);
s0M = kep2car(kep_Mc,mu_sun);

[tv,visibility_sun,visibility_mars] = evw_sun(s0E, state_ellisse(1,:), Tellisse*0.5, s0M);

figure()
hold on;
grid on;
plot(tv/3600,visibility_sun,'y','Linewidth',3);
plot(tv/3600,visibility_mars,'r','Linewidth',3);
xlabel('time of flight [hours]','Interpreter','latex');
ylabel('Earth link availability [1 = yes, 0 = no]','Interpreter','latex');
title('\textbf{Man 1 escape link availability}','Interpreter','latex');
legend('Location','southeast','Sun','Mars');

% Man 2 esc
kep_Ec = uplanet(date_escape + (Tellisse*0.5)/24/3600,3);
s0E = kep2car(kep_Ec,mu_sun);

kep_Mc = uplanet(date_escape + (Tellisse*0.5)/24/3600,4);
s0M = kep2car(kep_Mc,mu_sun);

[tv,visibility_sun,visibility_mars] = evw_sun(s0E, iperbole(1,:), Tellisse*5, s0M);

figure()
hold on;
grid on;
plot(tv/3600,visibility_sun,'y','Linewidth',3);
plot(tv/3600,visibility_mars,'r','Linewidth',3);
xlabel('time of flight [hours]','Interpreter','latex');
ylabel('Earth link availability [1 = yes, 0 = no]','Interpreter','latex');
title('\textbf{Man 2 escape link availability}','Interpreter','latex');
legend('Location','southeast','Sun','Mars');