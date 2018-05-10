clear 
close all
clc
%% DATA SECTION------------------------------------------------------------
addpath(genpath('./functions'))
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
medate_depi = [2028 09 1 12 00 00];      % Initial departure date
medate_depf = [2035 06 19 12 00 00];      % Final departure date
medate_arri = [2029 03 06 12 00 00];      % Initial arrival date
medate_arrf = [2035 05 10 12 00 00];      % Final arrival date

% Conversion to mjd2000----------------------------------------------------
medate_depi_mjd2000 = date2mjd2000(medate_depi);
medate_depf_mjd2000 = date2mjd2000(medate_depf);
medate_arri_mjd2000 = date2mjd2000(medate_arri);
medate_arrf_mjd2000 = date2mjd2000(medate_arrf);

% Definition of launch and arrival window arrays---------------------------
N = 5;                                              % Days per step 
medep_dates = medate_depi_mjd2000:N:medate_depf_mjd2000;
mearr_dates = medate_arri_mjd2000:N:medate_arrf_mjd2000;

%%
[mePCdata] = porkChopDatav3(planet_2,planet_1,medep_dates,mearr_dates);
dv1_dMAX=3;
porkChopPlot(mePCdata,1,0,1,dv1_dMAX);

%% Optimization throught ga------------------------------------------------
A = []; b = []; Aeq = []; beq = [];
lb = [medate_depi_mjd2000 medate_arri_mjd2000];
ub = [medate_depf_mjd2000 medate_arrf_mjd2000];
% options = optimoptions('ga','PopulationSize',100,...
%                             'MaxGenerations',20,...
%                             'FunctionTolerance',0,...
%                             'UseParallel',true,...
%                             'ConstraintTolerance',1e-10);
% sol_ga = ga(@(X)ffdv(X,2,planet_2,planet_1),2,A,b,Aeq,beq,lb,ub,@(X)dv1conme(X),options);

%% Optimization throught fmincon-------------------------------------------
%x0= [sol_ga(1) sol_ga(2)];
x0 = [10496 10832;
     11261 11592;
     12101 12332];
options_fmincon = optimoptions('fmincon','UseParallel',true);

dates=[];
DV=[];
for i = 1:size(x0,1)
dates(i,:) = fmincon(@(X)ffdv(X,2,planet_2,planet_1),x0(i,:),...
                      A,b,Aeq,beq,lb,ub,@(X)dv1conme(X),options_fmincon);
[dv_d,dv_a] = evaldv (dates(i,:),planet_2,planet_1);
DV(i,:) = [dv_d dv_a];
end

dv = [];
for i = 1:3
a = -mu_mars/DV(i,1)^2;
Vph = sqrt(2*mu_mars/(mr_mars+600)-mu_mars/a);
dv(i) = Vph-3.87;
end

