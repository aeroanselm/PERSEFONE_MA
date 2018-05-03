clear 
close all
clc
%% DATA SECTION------------------------------------------------------------
addpath(genpath('/home/alessandro/Documenti/GitHub/PERSEFONE_MA/functions'))
planet_1 = 3;                           % Earth number according to uplanet function
planet_2 = 4;                           % Mars number according to uplanet function
mu_earth = astroConstants(13);          % Earth's planetary constant 
mu_mars = astroConstants(14);           % Mars' planetary constant
mu_sun = astroConstants(4);             % Sun planetary constant
G = astroConstants(1);                  % Universal gravitational constant
AU = astroConstants(2);                 % Astronomic Unit
m_mars = mu_mars/G;                     % Mass of Mars
m_sun = mu_sun/G;                       % Mass of the Sun
m_earth = mu_earth/G;                   % Mass of the Earth
%m_phobos = 
mr_earth = astroConstants(23);          % Earth's mean radius
mr_mars = astroConstants(24);           % Mars' mean radius

%% Definition of the launch and arrival window----------------------------- 
date_depi = [2023 07 20 12 00 00];      % Initial departure date
date_depf = [2025 06 19 12 00 00];      % Final departure date
date_arri = [2024 11 06 12 00 00];      % Initial arrival date
date_arrf = [2026 05 10 12 00 00];      % Final arrival date

% Conversion to mjd2000----------------------------------------------------
date_depi_mjd2000 = date2mjd2000(date_depi);
date_depf_mjd2000 = date2mjd2000(date_depf);
date_arri_mjd2000 = date2mjd2000(date_arri);
date_arrf_mjd2000 = date2mjd2000(date_arrf);

% Definition of launch and arrival window arrays---------------------------
N = 5;                                              % Days per step 
dep_dates = date_depi_mjd2000:N:date_depf_mjd2000;
arr_dates = date_arri_mjd2000:N:date_arrf_mjd2000;

%% Optimization throught ga------------------------------------------------
A = []; b = []; Aeq = []; beq = [];
lb = [8601 9301];
ub = [9076 9626];
options = optimoptions('ga','PopulationSize',10,...
                            'MaxGenerations',10,...
                            'FunctionTolerance',0,...
                            'UseParallel',true,...
                            'ConstraintTolerance',1e-6);
sol_ga = ga(@(X)ffdv(X,2,planet_1,planet_2),2,A,b,Aeq,beq,lb,ub,@(X)dv1con(X),options);

%% Optimization throught fmincon-------------------------------------------
x0 = [sol_ga(1) sol_ga(2)];
options_fmincon = optimoptions('fmincon','UseParallel',true);
dates = fmincon(@(X)ffdv(X,2,planet_1,planet_2),x0,...
                      A,b,Aeq,beq,lb,ub,@(X)dv1con(X),options_fmincon);
[dv_d,dv_a] = evaldv (dates,planet_1,planet_2);

%% Transfer orbit definition------------------------------------------------
DAYd = dates(1,1);
DAYa = dates(1,2);
tof_sol=(DAYa-DAYd)*24*3600;
options=odeset('Reltol',1e-13,'AbsTol',1e-14);

% Earth
kep_Ed=uplanet(DAYd,planet_1);
state_Ed=kep2car(kep_Ed,mu_sun);
T_E=period(kep_Ed,mu_sun);
[tv_E,state_E]=ode113(@(t,y)dyn_2BP(t,y,mu_sun),[0,T_E],state_Ed,options);
kep_Ea=kep_Ed;
kep_Ea(6)=kep_Ed(6)+dt2dtheta(kep_Ed,mu_sun,tof_sol);
state_Ea=kep2car(kep_Ea,mu_sun);

% Mars
kep_Md=uplanet(DAYd,planet_2);
state_Md=kep2car(kep_Md,mu_sun);
T_M=period(kep_Md,mu_sun);
[tv_M,state_M]=ode113(@(t,y)dyn_2BP(t,y,mu_sun),[0,T_M],state_Md,options);
kep_Ma=kep_Md;
kep_Ma(6)=kep_Md(6)+dt2dtheta(kep_Md,mu_sun,tof_sol);
state_Ma=kep2car(kep_Ma,mu_sun);

[At,Pt,Et,~,V1] = lambertMR(state_Ed(1:3),state_Ma(1:3),tof_sol,mu_sun,0,0,0,0);

% Transfer
state0 = [state_Ed(1:3) V1];
kep_T=car2kep(state0,mu_sun);
T_T=period(kep_T,mu_sun);
options=odeset('Reltol',1e-13,'AbsTol',1e-14);
[t_tran,state_t]=ode113(@(t,y)dyn_2BP(t,y,mu_sun),[0,T_T],state0,options);

%% Capture hyperbola
mars_rsoi = norm(state_Ma(1:3))*(m_mars/m_sun)^(2/5);
[T_mars_t,state_Mt]=ode113(@(t,y)dyn_2BP(t,y,mu_sun),0:60:tof_sol,state_Md,options);
[T_t,state_tt]=ode113(@(t,y)dyn_2BP(t,y,mu_sun),T_mars_t,state0,options);
[rr_capture,vv_capture,date_capture] = patching(state_Mt,state_tt,T_t,DAYd,mars_rsoi);
ht = cross(rr_capture,vv_capture)/norm(cross(rr_capture,vv_capture));
S = vv_capture/norm(vv_capture);
B = cross(S,ht);
k = [0 0 1];                            % Unit vector normal to the ecliptic plane
T = cross(S,k)/norm(cross(S,k));
R = cross(S,T);
alfa = asin(dot(B,R));


