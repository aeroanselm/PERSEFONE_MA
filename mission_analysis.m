% mission_analysis.m
%  
% TYPE:
%   Main Script
%
% DESCRIPTION: 
%   This routine is used to obtain:
%       1 - DV budget for Mars transfer
%       2 - Minimum DV_tot
%       3 - Minimum constrained DV_tot
%       4 - Feasible trajectories
%       5 - Simil-Hohmann transfer 
%       6 - Optimization with genetic algorithm (ga)
%       7 - Optimization: ga -> fmincon
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
%   Creation: 15/04/2018
%
% CHANGELOG: 
%   ALESSANDRO MARIA MASSERINI creation of this script - 15/04/2018
%
% --------------------------------------------------------------------------

clear 
close all
clc
%% DATA SECTION
planet_1 = 3;                           % Earth number according to uplanet function
planet_2 = 4;                           % Mars number according to uplanet function
mu_earth = astroConstants(13);          % Earth's planetary constant 
mu_mars = astroConstants(14);           % Mars' planetary constant
mu_sun = astroConstants(4);             % Sun planetary constant
G = astroConstants(1);                  % Universal gravitational constant
m_mars = mu_mars/G;                     % Mass of Mars
m_sun = mu_sun/G;                       % Mass of the Sun
m_earth = mu_earth/G;                   % Mass of the Earth
mr_earth = astroConstants(23);          % Earth's mean radius
mr_mars = astroConstants(24);           % Mars' mean radius

%% Definition of the launch and arrival window 
date_depi = [2023 01 01 12 00 00];      % Initial departure date
date_depf = [2027 12 31 12 00 00];      % Final departure date
date_arri = [2024 03 01 12 00 00];      % Initial arrival date
date_arrf = [2030 06 31 12 00 00];      % Final arrival date

% Conversion to mjd200
date_depi_mjd2000 = date2mjd2000(date_depi);
date_depf_mjd2000 = date2mjd2000(date_depf);
date_arri_mjd2000 = date2mjd2000(date_arri);
date_arrf_mjd2000 = date2mjd2000(date_arrf);

% Definition of launch and arrival window arrays
N = 5;                                              % Days per step 
dep_dates = date_depi_mjd2000:N:date_depf_mjd2000;
arr_dates = date_arri_mjd2000:N:date_arrf_mjd2000;

% Martian seasons expressed as true anamalies
thetai_w = (-17-7.5)*pi/180;
thetaf_w = (-17+7.5)*pi/180;
thetac = [thetai_w thetaf_w];

% Pork chop plot of the Lambert's problem
[PCdata] = porkChopDatav3(planet_1,planet_2,dep_dates,arr_dates);
save PCdata.mat
dep=PCdata.dep;
arr=PCdata.arr;
TOF=PCdata.TOF;
dv1=PCdata.dv1;
dv2=PCdata.dv2;
DV=PCdata.DV;
dv1_dMAX=3.5;
%PCplot(dep,arr,TOF,DV,1,'Earth to Mars \DeltaV_{tot}','[Km/s]',10);
porkChopPlot(PCdata,1,0,2,dv1_dMAX);

%% Optimization throught fmincon
A = []; b = []; Aeq = []; beq = [];
lb = [date_depi_mjd2000 date_arri_mjd2000];
ub = [date_depf_mjd2000 date_arrf_mjd2000];
x0 = [9041 9356; 
      9801 10136];
flag_ff = 2;  
options = optimoptions('fmincon','UseParallel',true);
[n,m]=size(x0);
dates=zeros(n,m);
dv_d=zeros(n,1);
dv_a=zeros(n,1);

for it=1:2
    dates(it,:) = fmincon(@(X)ffdv(X,flag_ff,planet_1,planet_2),x0(it,:),A,b,Aeq,beq,lb,ub,@(X)dv1con(X),options);
    [dv_d(it),dv_a(it)] = evaldv (dates(it,:),planet_1,planet_2);
end

%% Transfer orbit definition
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

state0 = [state_Ed(1:3) V1];
kep_T=car2kep(state0,mu_sun);
T_T=period(kep_T,mu_sun);

options=odeset('Reltol',1e-13,'AbsTol',1e-14);
[t_tran,state_t]=ode113(@(t,y)dyn_2BP(t,y,mu_sun),[0,T_T],state0,options);

% Plot
figure();
hold on
%grid on
axis equal
set(gca,'Color','k')
drawPlanet('Sun',[0 0 0],gca,20);
drawPlanet('Earth',state_Ed(1:3),gca,1000);
drawPlanet('Mars',state_Ma(1:3),gca,1000);
p1=plot3(state_E(:,1),state_E(:,2),state_E(:,3),'--b');
p2=plot3(state_M(:,1),state_M(:,2),state_M(:,3),'--w');
p3=plot3(state_t(:,1),state_t(:,2),state_t(:,3),'--r','LineWidth',2);
%p4=plot3(state_Tp(:,1),state_Tp(:,2),state_Tp(:,3),'r','LineWidth',2);
p5=plot3(state_Ed(1),state_Ed(2),state_Ed(3),'*w');
p6=plot3(state_Ea(1),state_Ea(2),state_Ea(3),'ow');
plot3(state_Md(1),state_Md(2),state_Md(3),'*w');
plot3(state_Ma(1),state_Ma(2),state_Ma(3),'ow');
%legend([p1 p2 p3 p5 p6],'Earth orbit','Mars orbit','Transfer orbit', 'Start','Finish');
title('Trajectories');
xlabel('x [Km]');
ylabel('y [Km]');
zlabel('z [Km]');

%% SAA - Sun Aspect Angle evolution
[theta,S] = saa(mu_sun,state0(1:3),state0(4:6),tof_sol);

figure()
hold on 
grid on
axis equal
plot(1:length(theta),theta)

%% Escape hyperbola
% earth_rsoi = norm(state_Ed(1:3))*(m_earth/m_sun)^(2/5); 
% [T_earth_t,state_Et]=ode113(@(t,y)dyn_2BP(t,y,mu_sun),0:1:4*24*3600,state_Ed,options);
% [T_t,state_t]=ode113(@(t,y)dyn_2BP(t,y,mu_sun),T_earth_t,state0,options);
% [rr_escape,vv_escape,date_escape] = patching(state_Et,state_t,T_t,DAYd,earth_rsoi);
% escape_time = date_escape-DAYd;         % Expressed in days

%% Capture hyperbola
mars_rsoi = norm(state_Ma(1:3))*(m_mars/m_sun)^(2/5);
[T_mars_t,state_Mt]=ode113(@(t,y)dyn_2BP(t,y,mu_sun),0:60:tof_sol,state_Md,options);
[T_t,state_t]=ode113(@(t,y)dyn_2BP(t,y,mu_sun),T_mars_t,state0,options);
[rr_capture,vv_capture,date_capture] = patching(state_Mt,state_t,T_t,DAYd,mars_rsoi);
ht = cross(rr_capture,vv_capture)/norm(cross(rr_capture,vv_capture));
S = vv_capture/norm(vv_capture);
B = cross(S,ht);
k = [0 0 1];                            % Unit vector normal to the ecliptic plane
T = cross(S,k)/norm(cross(S,k));
R = cross(S,T);
alfa = asin(dot(B,R));

%% Earth visibility window from 2025/08/31 08:00 to 2026/03/31 08:00
load marseph20250831_20260331.mat
load phoboseph20250831_20260331.mat
S_Mars = marseph20250831_20260331;
S_Phobos = phoboseph20250831_20260331;
dt = date2mjd2000([2025 08 31 08 05 00])-date2mjd2000([2025 08 31 08 00 00]); % 5minutes
dates = [date2mjd2000([2025 08 31 08 00 00]):dt:date2mjd2000([2025 08 31 08 05 00])];
vd = evw(mr_mars, S_Mars, S_Phobos,dates,1,100);




























% figure()
% hold on
% grid on
% axis equal
% drawPlanet('Sun',[0 0 0],gca,20);
% drawPlanet('Earth',state_Ed(1:3),gca,1000);
% drawPlanet('Mars',state_Ma(1:3),gca,100);
% drawPlanet('Earth',state_Ea(1:3),gca,100);
% drawPlanet('Mars',state_Md(1:3),gca,1000);
% plot3(state_t(:,1), state_t(:,2), state_t(:,3),'--');

%% Plot trajectories











%% Optimization throught genetic algorithm
































































% 
% 
% 
% 
% 
% 
% kep_E = uplanet(date_depi_mjd2000,planet_1);
% %kep_E(end) = 0;
% %kep_E(6)*180/pi
% kep_M = uplanet(date_depi_mjd2000,planet_2);
% kep_M(6)*180/pi
% 
% state_E = kep2car(kep_E,mu_sun);
% state_M = kep2car(kep_M,mu_sun);
% 
% r_E = state_E(1:3);
% r_M = state_M(1:3);
% 
% %alpha = acos(dot(r_E,r_M)/(norm(state_E)*norm(state_M)))*180/pi
% 
% ratio = 1500;
% 
% figure()
% hold on
% grid on
% axis equal
% drawPlanet('Sun',[0 0 0],gca,20);
% drawPlanet('Earth',r_E,gca,ratio);
% drawPlanet('Mars',r_M,gca,ratio);
% plot(norm(r_E),0,'x')
% 
% 
