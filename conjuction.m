clear 
close all
clc
%% DATA SECTION------------------------------------------------------------
addpath(genpath('./functions'))         % Adding the functions folder to path [linux version, works even on Windows]
planet_1 = 3;                           % Earth number according to uplanet function
planet_2 = 4;                           % Mars number according to uplanet function
mu_earth = astroConstants(13);          % Earth's planetary constant 
mu_mars = astroConstants(14);           % Mars' planetary constant
mu_sun = astroConstants(4);             % Sun planetary constant
mu_phobos = 0.0007112;                  % Phobos gravitational constant
G = astroConstants(1);                  % Universal gravitational constant
AU = astroConstants(2);                 % Astronomic Unit
m_mars = mu_mars/G;                     % Mass of Mars
m_sun = mu_sun/G;                       % Mass of the Sun
m_earth = mu_earth/G;                   % Mass of the Earth
m_phobos = mu_phobos/G;                 % Mass of Phobos
mr_earth = astroConstants(23);          % Earth's mean radius
mr_mars = astroConstants(24);           % Mars' mean radius
options = odeset('Reltol',1e-13,'AbsTol',1e-14);

day_arrival = [2025 08 31 08 32 13];
day_arrival = date2mjd2000(day_arrival);

kep_Ea = uplanet(day_arrival, planet_1);
kep_Ma = uplanet(day_arrival, planet_2);

S_Ea = kep2car(kep_Ea, mu_sun);
S_Ma = kep2car(kep_Ma, mu_sun);

%% Propagation of motion
day = 3600*24;
t = day;
state_E = S_Ea;
state_M = S_Ma;
n_days = 0;

while(checkfcn(state_E,state_M))
    n_days = n_days + 1;
%     [~,state_E] = ode113(@(t,x)dyn_2BP(t,x,mu_sun),[0 t],state_E,options);
%     [~,state_M] = ode113(@(t,x)dyn_2BP(t,x,mu_sun),[0 t],state_M,options);
    kep_E = uplanet(day_arrival + n_days, planet_1);
    kep_M = uplanet(day_arrival + n_days, planet_2);
    state_E = kep2car(kep_E, mu_sun);
    state_M = kep2car(kep_M, mu_sun);
    state_E = state_E(end,:);
    state_M = state_M(end,:);
    
end

figure()
hold on;
grid on;
axis equal;
title('\textbf{Mars Conjunction}','Interpreter','latex');
xlabel('\textbf{[km]}','Interpreter','latex');
ylabel('\textbf{[km]}','Interpreter','latex');
zlabel('\textbf{[km]}','Interpreter','latex');

p1 = plot3(0, 0, 0,'*y');
p2 = plot3(S_Ea(1), S_Ea(2), S_Ea(3),'*b');
p3 = plot3(S_Ma(1), S_Ma(2), S_Ma(3),'*r');
p4 = plot3(state_E(1), state_E(2), state_E(3),'ob');
p5 = plot3(state_M(1), state_M(2), state_M(3),'or');

n_days_conj = 0;

while(checkfcn(state_E,state_M) == false)
    n_days_conj = n_days_conj + 1;
%     [~,state_E] = ode113(@(t,x)dyn_2BP(t,x,mu_sun),[0 t],state_E,options);
%     [~,state_M] = ode113(@(t,x)dyn_2BP(t,x,mu_sun),[0 t],state_M,options);
    kep_E = uplanet(day_arrival + n_days + n_days_conj, planet_1);
    kep_M = uplanet(day_arrival + n_days + n_days_conj, planet_2);
    state_E = kep2car(kep_E, mu_sun);
    state_M = kep2car(kep_M, mu_sun);
    state_E = state_E(end,:);
    state_M = state_M(end,:);
end

p6 = plot3(state_E(1), state_E(2), state_E(3),'xb');
p7 = plot3(state_M(1), state_M(2), state_M(3),'xr');
xlim([-3e+8 3e+8]);
ylim([-3e+8 3e+8]);
zlim([-3e+8 3e+8]);
legend([p1 p2 p3 p4 p5 p6 p7],'Location','bestoutside',...
    'Sun','Arrival (Earth)', 'Arrival (Mars)','Start of conjunction (Earth)',...
    'Start of conjunction (Mars)', 'End of conjunction (Earth)',...
    'End of conjunction (Mars)');

day_conj_start = day_arrival + n_days;
day_conj_end = day_conj_start + n_days_conj;

date_conj_start = mjd20002date(day_conj_start);
date_conj_end = mjd20002date(day_conj_end);

str_start = sprintf('%d ', date_conj_start(1:3));
fprintf('Conjunction starts on: %s\n', str_start);
str_end = sprintf('%d ', date_conj_end(1:3));
fprintf("Conjunction ends on: %s\n", str_end)

function [pippo] = checkfcn (state_E, state_M)
    SE = norm(state_E(1:3));
    SM = norm(state_M(1:3));
    EM = norm(state_M(1:3) - state_E(1:3));
    
    p = (SE + SM + EM)/2;
    A = sqrt(p*(p - SE)*(p - SM)*(p - EM));
    SH = 2*A/EM;
    
    if SH <= SM*sin(1.7*pi/180)
        pippo = false;
    else
        pippo = true;
        return 
    end
end


