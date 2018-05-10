function [] = transforb(dates,planet_1, planet_2, mu_sun)

DAYd = dates(1,1);          % Selected departure date [mjd2000]
DAYa = dates(1,2);          % Selected arrival date   [mjd2000]
tof_sol = (DAYa - DAYd)*24*3600;        % Transfer time

% Setting ode options
options = odeset('Reltol',1e-13,'AbsTol',1e-14);

% Earth - State definition at departure and arrival dates
kep_Ed = uplanet(DAYd,planet_1);
state_Ed = kep2car(kep_Ed,mu_sun);
T_E = period(kep_Ed,mu_sun);
[tv_E,state_E] = ode113(@(t,y)dyn_2BP(t,y,mu_sun),[0,T_E],state_Ed,options);
kep_Ea = kep_Ed;
kep_Ea(6) = kep_Ed(6)+dt2dtheta(kep_Ed,mu_sun,tof_sol);
state_Ea = kep2car(kep_Ea,mu_sun);

% Mars - State definition at departure and arrival dates
kep_Md = uplanet(DAYd,planet_2);
state_Md = kep2car(kep_Md,mu_sun);
T_M = period(kep_Md,mu_sun);
[tv_M,state_M] = ode113(@(t,y)dyn_2BP(t,y,mu_sun),[0,T_M],state_Md,options);
kep_Ma = kep_Md;
kep_Ma(6) = kep_Md(6)+dt2dtheta(kep_Md,mu_sun,tof_sol);
state_Ma = kep2car(kep_Ma,mu_sun);

% Transfer arc computation
[~,~,~,~,V1,~] = lambertMR(state_Ed(1:3),state_Ma(1:3),tof_sol,mu_sun,0,0,0,0);

% Transfer orbit computation
state0 = [state_Ed(1:3) V1];        % Initial state vector
kep_T = car2kep(state0,mu_sun);     % Keplerian elements of the transfer orbit
T_T = period(kep_T,mu_sun);         % Transfer orbit period

[tv_T,state_T] = ode113(@(t,y)dyn_2BP(t,y,mu_sun),[0,T_T],state0,options);
