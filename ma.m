% mission_analysisV2.m
%  
% TYPE:
%   Main Script
%
% DESCRIPTION: 
%   This routine is used to obtain:
%       1 - DV budget for Mars system insertion
%       2 - Mamoeuvres optimization to approach the orbit of Phobos
%       3 - Mars-Earth leg design
%       4 - Sun Aspect Angle evaluation for each significant phase 
%       5 - Earth visibility windows evaluation for each significant phase 
%       6 - Evaluation of the beta angle for TCS support
%       END - Final plot and considerations
%
%   Note:
%
% NOTATION:
%
% REFERENCES:
%   - Author,"Book's Name or Article", dd/mm/yyyy.
%
% FUTURE DEVELOPMENT:
%   Work in progress.
%
% ORIGINAL VERSION:
%   17/12/2016, ALESSANDRO MARIA MASSERINI
%
% AUTHORS:
%
%   Name: ALESSANDRO MARIA 
%   Surname: MASSERINI
%   ID number: 808451
%   Contact: alessandro.masserini@mail.polimi.it
%
%   Course: Space Engineering
%   Department: DAER
%   University: Politecnico di Milano
%   Class: Space Mission Analysis And Design 
%   Creation: 11/05/2018
%
% CHANGELOG: 
%   ALESSANDRO MARIA MASSERINI creation of this script - 11/05/2018
%
% -------------------------------------------------------------------------
clear 
close all
clc

%% DATA SECTION------------------------------------------------------------
addpath(genpath('./functions'))         % Adding the functions folder to path [linux version, works even on Windows]
ma_data;

%% Definition of the launch and arrival window----------------------------- 
date_depi = [2023 07 01 12 00 00];      % Initial departure date
date_depf = [2025 09 30 12 00 00];      % Final departure date
date_arri = [2024 05 01 12 00 00];      % Initial arrival date
date_arrf = [2026 08 31 12 00 00];      % Final arrival date

% Conversion to mjd2000----------------------------------------------------
date_depi_mjd2000 = date2mjd2000(date_depi);    % Initial departure date [mjd2000]
date_depf_mjd2000 = date2mjd2000(date_depf);    % Final departure date   [mjd2000]
date_arri_mjd2000 = date2mjd2000(date_arri);    % Initial arrival date   [mjd2000]
date_arrf_mjd2000 = date2mjd2000(date_arrf);    % Final arrival date     [mjd2000]

% Definition of launch and arrival window arrays---------------------------
N = 5;                                              % Days per step 
dep_dates = date_depi_mjd2000:N:date_depf_mjd2000;  % Departure dates vector
arr_dates = date_arri_mjd2000:N:date_arrf_mjd2000;  % Arrival dates vector

%% Pork chop plot of the Lambert's problem----------------------------------
[PCdata] = porkChopDatav3(planet_1,planet_2,dep_dates,arr_dates);
save PCdata.mat
dep=PCdata.dep;
arr=PCdata.arr;
TOF=PCdata.TOF;
dv1=PCdata.dv1;
dv2=PCdata.dv2;
DV=PCdata.DV;
dv1_dMAX=3.5;
porkChopPlot(PCdata,1,0,1,dv1_dMAX);

%% Optimization throught fmincon-------------------------------------------
A = []; b = []; Aeq = []; beq = [];
lb = [date_depi_mjd2000 date_arri_mjd2000];     % Lower boundary f the ga domain
ub = [date_depf_mjd2000 date_arrf_mjd2000];     % Upper boundary f the ga domain
x0 = [9040.38764537296 9373.85570383965]; 

% Setting fmincon options 
options_fmincon = optimoptions('fmincon','UseParallel',true);

% Calling 'fmincon'
dates = fmincon(@(X)ffdv(X,2,planet_1,planet_2),x0,...
                      A,b,Aeq,beq,lb,ub,@(X)dv1con(X),options_fmincon);
[dv_d,dv_a] = evaldv (dates,planet_1,planet_2); % Evaluation of the fmincon solution

%% Keplerian elements of Phobos Orbit at 2025/08/31
kep_p = [9.379005712581698E+03 1.547082708359828E-02 (2.778424215200678E+01)*pi/180 ...
        (8.299695327985185E+01)*pi/180 (1.052549189428009E+02)*pi/180 (2.122402357425826E+01)*pi/180];
		
kep_p = [kep_p(1) 0 kep_p(3) kep_p(4) kep_p(5) kep_p(6)];

%% Transfer orbit definition------------------------------------------------
DAYd = dates(1);
DAYa = dates(2);
save ('./data/departure.mat','DAYd');
save ('./data/arrival.mat','DAYa');
[state_T, state_1, state_2, tof_sol] = transforb(dates,planet_1, planet_2, mu_sun);
kepseason = uplanet(DAYa,4);
[season,Ls] = mseasons(kepseason(6)*180/pi);
%% Capture hyperbola
state_Md = state_2(1,:);
state_Ma = state_2(end,:);
state0 = state_T(1,:);
mars_rsoi = norm(state_Ma(1:3))*(m_mars/m_sun)^(2/5);
[tv_Mt,state_Mt]=ode113(@(t,y)dyn_2BP(t,y,mu_sun),0:60:tof_sol,state_Md,options);
[tv_t,state_tt]=ode113(@(t,y)dyn_2BP(t,y,mu_sun),tv_Mt,state0,options);
[rr_capture,vv_capture,date_capture] = patching(state_Mt,state_tt,tv_t,DAYd,mars_rsoi);
save('./data/DAYc.mat','date_capture');
% B-plane definition
ht = cross(rr_capture,vv_capture)/norm(cross(rr_capture,vv_capture));   % Angular momentum unit vector of the transfer orbit
[R,T,S] = bplane(vv_capture);
%B = cross(S,ht);                                                        % Unit vector of the impact parameter

%% Impact parameter selection
theta_deg = 30;
theta = theta_deg*pi/180;
delta =  7924;
B = sin(theta)*R + cos(theta)*T;
Delta = delta*B;    
gamma = asin(norm(Delta)/mars_rsoi);
rrpar = mars_rsoi*cos(gamma)*-S;
rr_capture = rrpar+Delta;
state_h = [rr_capture vv_capture];
rp_f = kep_p(1);
[dvv, dv, dt,stateman,rp_h, kepf, sol] = mancombo(state_h,mu_mars,kep_p(3),kep_p(4),0.9,rp_f)
s1 = stateman.state_1;
s2 = stateman.state_2;
s3 = stateman.state_3;
s4 = stateman.state_4;
s5 = stateman.state_5;
h_closest = rp_h - mr_mars;
save ('Phobos_Orbit.mat','kepf');
save('./data/iperbole.mat','s1');
save('./data/ellisse.mat','s2');
save('./data/cps.mat','s3');
save('./data/cOM.mat','s4');
save('./data/circ.mat','s5');
save('./data/dt.mat','dt');
%% Only for debug
delta_range = [7924];%, 7924, 8924];
theta_range_deg = 20:80;
theta_range = theta_range_deg*pi/180;

str = [];
figure()
hold on;
grid on;
title('\DeltaV vs \theta');
for i = 1:length(delta_range)
	delta = delta_range(i);
	dvtot = [];
	for j = 1:length(theta_range)
		B = sin(theta_range(j))*R + cos(theta_range(j))*T;
		Delta = delta*B;
		gamma = asin(norm(Delta)/mars_rsoi);
		rrpar = mars_rsoi*cos(gamma)*-S;
		rr_capture = rrpar+Delta;
		state_h = [rr_capture vv_capture];
		rp_f = kep_p(1);
		[dvv, dv, dt,stateman,rp_h, kepf, sol] = mancombo(state_h,mu_mars,kep_p(3),kep_p(4),0.9,rp_f);
		dvtot(j) = sum(dv);
    end
	plot(theta_range*180/pi,dvtot);
    % insert legend
end

%% Plot
figure()
hold on
axis equal
title('\textbf{Mars system insertion}','Interpreter','latex');
drawPlanet('Mars',[0 0 0],gca,1);
p1 = plot3(s1(:,1),s1(:,2),s1(:,3),'Linewidth',1.5);
p2 = plot3(s2(:,1),s2(:,2),s2(:,3),'Linewidth',1.5);
p3 = plot3(s3(:,1),s3(:,2),s3(:,3),'Linewidth',1.5);
p4 = plot3(s4(:,1),s4(:,2),s4(:,3),'Linewidth',1.5);
p5 = plot3(s5(:,1),s5(:,2),s5(:,3),'Linewidth',1.5);
legend([p1 p2 p3 p4 p5],'Location','southeast','Hyperbolic leg',...
        'HEO','i and a change',...
        '\Omega change','Circularize');

%% Lambert arc recalculation
kep_Msoi = uplanet(date_capture,4);
state_Msoi = kep2car(kep_Msoi, mu_sun);
state_Ed = state_1(1,:);
[~,~,~,~,~,V2] = lambertMR(state_Ed(1:3),rr_capture+state_Msoi(1:3),tof_sol-(DAYa-date_capture)*3600*24,mu_sun,0,0,0,0);
V_relsoi=V2-state_Msoi(4:6);
dsm_a = norm(vv_capture-V_relsoi);

%% Return to Earth
% Definition of the launch and arrival window----------------------------- 
medate_depi = [2024 01 01 12 00 00];      % Initial departure date
medate_depf = [2035 12 31 12 00 00];      % Final departure date
medate_arri = [2024 01 01 12 00 00];      % Initial arrival date
medate_arrf = [2035 12 31 12 00 00];      % Final arrival date

% Conversion to mjd2000----------------------------------------------------
medate_depi_mjd2000 = date2mjd2000(medate_depi);
medate_depf_mjd2000 = date2mjd2000(medate_depf);
medate_arri_mjd2000 = date2mjd2000(medate_arri);
medate_arrf_mjd2000 = date2mjd2000(medate_arrf);

% Definition of launch and arrival window arrays---------------------------
N = 5;                                              % Days per step 
medep_dates = medate_depi_mjd2000:N:medate_depf_mjd2000;
mearr_dates = medate_arri_mjd2000:N:medate_arrf_mjd2000;
% medep_dates = 10000:N:11830;
% mearr_dates = 10790:N:12020;
%%
[mePCdata] = porkChopDatav3(planet_2,planet_1,medep_dates,mearr_dates);
dv1_dMAX=3;
porkChopPlot(mePCdata,1,0,1,dv1_dMAX);

%% Optimization throught fmincon-------------------------------------------
lb = [medate_depi_mjd2000 medate_arri_mjd2000];
ub = [medate_depf_mjd2000 medate_arrf_mjd2000];
%x0= [sol_ga(1) sol_ga(2)];
% x0 = [10496 10832;
%      11261 11592;
%      12101 12332];

x0 = [11261 11592];
options_fmincon = optimoptions('fmincon','UseParallel',true);


dates = fmincon(@(X)ffdv(X,2,planet_2,planet_1),x0,...
                      A,b,Aeq,beq,lb,ub,@(X)dv1conme(X),options_fmincon);
[dv_d,dv_a] = evaldv (dates,planet_2,planet_1);
DV = [dv_d dv_a];
a = -mu_mars/DV(1)^2;
Vph = sqrt(2*mu_mars/(mr_mars+600)-mu_mars/a);
dv= Vph-3.87;

%% Transfer orbit definition
meDAYd = dates(1);
meDAYa = dates(2);
medates = [meDAYd, meDAYa];
save('./data/medates.mat','medates');
[mestate_T, mestate_M, mestate_E, metof_sol] = transforb(medates,planet_2, planet_1, mu_sun);


%% Escape hyperbola
mestate_Md = mestate_M(1,:);
mestate0 = mestate_T(1,:);
[metv_Mt,mestate_Mt]=ode113(@(t,y)dyn_2BP(t,y,mu_sun),0:60:metof_sol,mestate_Md,options);
[metv_t,mestate_tt]=ode113(@(t,y)dyn_2BP(t,y,mu_sun),metv_Mt,mestate0,options);
[rr_escape,vv_escape,date_escape] = patching(mestate_Mt,mestate_tt,metv_t,meDAYd,mars_rsoi);
save('./data/dateesc.mat','date_escape');
% B-plane definition
ht = cross(rr_escape,vv_escape)/norm(cross(rr_escape,vv_escape));   % Angular momentum unit vector of the transfer orbit
[R,T,S] = bplane(vv_escape);

%B = cross(S,ht);                                                        % Unit vector of the impact parameter
a_esc = -mu_mars/(norm(vv_escape)^2-2*mu_mars/mars_rsoi);
rp_esc = mr_mars+600;
vp_esc = sqrt(2*mu_mars/rp_esc-mu_mars/a_esc);
vinf_esc = sqrt(-mu_mars/a_esc);
delta_esc = mu_mars/vinf_esc^2*((1+vinf_esc^2*rp_esc/mu_mars)^2 - 1)^0.5;

theta_escdeg = 25;
theta_esc = theta_escdeg*pi/180;
B = sin(theta_esc)*R + cos(theta_esc)*T;
Delta_esc = delta_esc*B;
gammaesc = asin(delta_esc/mars_rsoi);
rrparesc = mars_rsoi*cos(gammaesc)*S;
rr_escape = Delta_esc + rrparesc;
state_esc = [rr_escape vv_escape];
kep_esc = car2kep(state_esc, mu_mars);
i_esc = kep_esc(3);
OM_esc = kep_esc(4);

[dv_cpesc, kep2esc] = cplane(kepf,OM_esc,i_esc,mu_mars);
dvesc = dv_cpesc;

%% LMO insertion and escape
v_cir = sqrt(mu_mars/kep2esc(1)); 
ra_esc = kep2esc(1);
a_eesc = (rp_esc + ra_esc)/2;
e_esc = (-rp_esc + ra_esc)/(rp_esc + ra_esc);
P_esc = a_eesc*(1 - e_esc^2);
Va = sqrt(mu_mars/P_esc)*(1-e_esc);
Vp = sqrt(mu_mars/P_esc)*(1+e_esc);
dv1h = Va - v_cir;
Vlmo = sqrt(mu_mars/rp_esc);
dv2h = Vlmo - Vp;

dv_e_h = vp_esc - Vp;

%% Plot escape
state_phobos = kep2car(kepf,mu_mars);
T_phobos = period(kepf,mu_mars);
[~,sphobos]=ode113(@(t,y)dyn_2BP(t,y,mu_mars),[0 T_phobos],state_phobos,options);
statecp = kep2car(kep2esc,mu_mars);
[~,scp]=ode113(@(t,y)dyn_2BP(t,y,mu_mars),[0 T_phobos],statecp,options);
kep_ellisse = [a_eesc e_esc kep2esc(3) kep_esc(4) kep_esc(5) pi];
Tellisse = period(kep_ellisse,mu_mars);
save('./data/Tellisseesc.mat','Tellisse');
state_ellisse = kep2car(kep_ellisse,mu_mars);
save('./data/ellissesc.mat','state_ellisse');
[~,ellisse]=ode113(@(t,y)dyn_2BP(t,y,mu_mars),[0 0.5*Tellisse],state_ellisse,options);
keplmo = kep_ellisse;
keplmo(1) = rp_esc; keplmo(2) = 0;
statelmo = kep2car(keplmo,mu_mars);
save('./data/sLMO.mat','statelmo');
Tlmo = period(keplmo,mu_mars);
save('./data/Tlmo','Tlmo');
[~,lmo]=ode113(@(t,y)dyn_2BP(t,y,mu_mars),[0 Tlmo],statelmo,options);
s_ip = ellisse(end,:);
vv_ip = vp_esc*(s_ip(4:6)/norm(s_ip(4:6)));
s_ip = [s_ip(1:3) vv_ip];
[~,iperbole]=ode113(@(t,y)dyn_2BP(t,y,mu_mars),[0 3600*4],s_ip,options);
save('./data/hypesc.mat','iperbole');
%%
figure()
hold on
axis equal
title('\textbf{Mars system escape}','Interpreter','latex');
drawPlanet('Mars',[0 0 0],gca,1);
q1 = plot3(sphobos(:,1),sphobos(:,2), sphobos(:,3),'Linewidth',1.5);
q2 = plot3(scp(:,1),scp(:,2), scp(:,3),'Linewidth',1.5);
q3 = plot3(ellisse(:,1),ellisse(:,2), ellisse(:,3),'Linewidth',1.5);
q4 = plot3(lmo(:,1),lmo(:,2), lmo(:,3),'Linewidth',1.5);
q5 = plot3(iperbole(:,1),iperbole(:,2), iperbole(:,3),'Linewidth',1.5);
legend([q1 q2 q3 q4 q5],'Location','southeast','Phobos orbit','i and \Omega change','Transfer ellipse','LMO','Hyberbolic escape');
%% Lambert arc recalculation
mekep_Msoi = uplanet(date_escape,4);
mestate_Msoi = kep2car(mekep_Msoi, mu_sun);
mekep_Ea = uplanet(meDAYa,3);
mestate_Ea = kep2car(mekep_Ea,mu_sun);
[~,~,~,~,V1,~] = lambertMR(mestate_Msoi(1:3)+rr_escape,mestate_Ea(1:3),(meDAYa-date_escape)*3600*24,mu_sun,0,0,0,0);
V_relsoi=V1-mestate_Msoi(4:6);
medsm_d = norm(vv_escape-V_relsoi);

% kep_esc(6) = 0;
% state_hp = kep2car(kep_esc,mu_mars);
% vvp_esc_u = state_hp(4:6)/norm(state_hp(4:6)); 
% VVp = Vp*vvp_esc_u;
% state_eesc = [state_hp(1:3) VVp];
% kep_eesc = car2kep(state_eesc, mu_mars);
% T_eesc = period(kep_eesc,mu_mars);
% [~,ellisse]=ode113(@(t,y)dyn_2BP(t,y,mu_mars),[0 T_eesc],state_eesc,options);
% state_lmo = [state_hp(1:3) Vlmo*vvp_esc_u];
% kep_lmo = car2kep(state_lmo,mu_mars);
% Tlmo = period(kep_lmo,mu_mars);
% [~,lmo]=ode113(@(t,y)dyn_2BP(t,y,mu_mars),[0 Tlmo],state_lmo,options);
% [~,iperbole]=ode113(@(t,y)dyn_2BP(t,y,mu_mars),[0 3600*8],state_hp,options);
