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
G = astroConstants(1);                  % Universal gravitational constant
AU = astroConstants(2);                 % Astronomic Unit
m_mars = mu_mars/G;                     % Mass of Mars
m_sun = mu_sun/G;                       % Mass of the Sun
m_earth = mu_earth/G;                   % Mass of the Earth
mu_phobos = 0.0007112;                  % Phobos gravitational constant
m_phobos = mu_phobos/G;                 % Mass of Phobos
mr_earth = astroConstants(23);          % Earth's mean radius
mr_mars = astroConstants(24);           % Mars' mean radius


%% Definition of the launch and arrival window----------------------------- 
date_depi = [2023 07 20 12 00 00];      % Initial departure date
date_depf = [2025 06 19 12 00 00];      % Final departure date
date_arri = [2024 11 06 12 00 00];      % Initial arrival date
date_arrf = [2026 05 10 12 00 00];      % Final arrival date

% Conversion to mjd2000----------------------------------------------------
date_depi_mjd2000 = date2mjd2000(date_depi);    % Initial departure date [mjd2000]
date_depf_mjd2000 = date2mjd2000(date_depf);    % Final departure date   [mjd2000]
date_arri_mjd2000 = date2mjd2000(date_arri);    % Initial arrival date   [mjd2000]
date_arrf_mjd2000 = date2mjd2000(date_arrf);    % Final arrival date     [mjd2000]

% Definition of launch and arrival window arrays---------------------------
N = 5;                                              % Days per step 
dep_dates = date_depi_mjd2000:N:date_depf_mjd2000;  % Departure dates vector
arr_dates = date_arri_mjd2000:N:date_arrf_mjd2000;  % Arrival dates vector

%% Optimization throught ga------------------------------------------------
A = []; b = []; Aeq = []; beq = [];
lb = [date_depi_mjd2000 date_arri_mjd2000];     % Lower boundary f the ga domain
ub = [date_depf_mjd2000 date_arrf_mjd2000];     % Upper boundary f the ga domain

% Setting ga options
% options = optimoptions('ga','PopulationSize',100,...
%                             'MaxGenerations',20,...
%                             'FunctionTolerance',0,...
%                             'UseParallel',true,...
%                             'ConstraintTolerance',1e-10);   

% Calling 'ga'                            
%sol_ga = ga(@(X)ffdv(X,2,planet_1,planet_2),2,A,b,Aeq,beq,lb,ub,@(X)dv1con(X),options);

%% Optimization throught fmincon-------------------------------------------
% x0 = [sol_ga(1) sol_ga(2)];                     % Initial guess is the ga solution
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
    
%% Transfer orbit definition------------------------------------------------
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

%% Capture hyperbola
mars_rsoi = norm(state_Ma(1:3))*(m_mars/m_sun)^(2/5);
[tv_Mt,state_Mt]=ode113(@(t,y)dyn_2BP(t,y,mu_sun),0:60:tof_sol,state_Md,options);
[tv_t,state_tt]=ode113(@(t,y)dyn_2BP(t,y,mu_sun),tv_Mt,state0,options);
[rr_capture,vv_capture,date_capture] = patching(state_Mt,state_tt,tv_t,DAYd,mars_rsoi);

ht = cross(rr_capture,vv_capture)/norm(cross(rr_capture,vv_capture));   % Angular momentum unit vector of the transfer orbit
S = vv_capture/norm(vv_capture);                                        % Velocity unit vector at the entry of the SOI
B = cross(S,ht);                                                        % Unit vector of the impact parameter
k = [0 0 1];                                                            % Unit vector normal to the ecliptic plane
T = cross(S,k)/norm(cross(S,k));                                        % Unit vector lying on B-plane and parallel wrt the ecliptic plane
R = cross(S,T);                                                         % Unit vector lying on B-plane and normal wrt S and T
theta0 = asin(dot(B,R));                                                % Angle between B and T 

% %% Selecting Ipact Parameter and Periapsis Radius to design the capture leg of the hyperbola
% theta_deg = 0:40;
% theta = theta_deg*pi/180;
% delta =  7924;
% [n,m] = size(delta);
% 
% figure()
% hold on
% grid on
% title('\Delta V depending on aiming radius vector');
% xlabel('Impact parameter [km]');
% ylabel('\Delta V_{tot} [km/s]');
% 
% dvtot = zeros(n,m);
% for k = 1:length(theta)
%     
%         B = sin(theta(k))*R + cos(theta(k))*T;
%         Delta = delta*B;    
%         gamma = asin(norm(Delta)/mars_rsoi);
%         rrpar = mars_rsoi*cos(gamma)*-S;
%         rr_capture = rrpar+Delta;
%         state_h = [rr_capture vv_capture];
%         rp_f = kep_p(1);
%         [dvv, dv, dt,stateman,rp_h, kepf, sol] = mancombo(state_h,mu_mars,kep_p(3),kep_p(4),0.9,rp_f);
%         
%         h_closest = rp_h - mr_mars;
%         if h_closest < 120
%              dvtot(k) = NaN;
%         else
%             dvtot(k) = sum(dv);
%         end
%         plot(theta,dvtot);
% %         str = sprintf("theta = %d", theta*180/pi);
% %         legend(str);
% end

%%
%plot(theta*180/pi,dvtot);

%%
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
[dvv, dv, dt,stateman,rp_h, kepf, sol] = mancombo(state_h,mu_mars,kep_p(3),kep_p(4),0.9,rp_f);
        
h_closest = rp_h - mr_mars;
%%
figure()
hold on
axis equal
drawPlanet('Mars',[0 0 0],gca,1);
comet3(stateman(:,1),stateman(:,2),stateman(:,3));%,'Linewidth',1.5)    
title('Approaching Phobos orbit manoeuvres')
legend('Mars','Trajectory')

%% Lambert arc recalculation
kep_Msoi = uplanet(date_capture,4);
state_Msoi = kep2car(kep_Msoi, mu_sun);

[~,~,~,~,~,V2] = lambertMR(state_Ed(1:3),rr_capture+state_Msoi(1:3),tof_sol-(DAYa-date_capture)*3600*24,mu_sun,0,0,0,0);
V_relsoi=V2-state_Msoi(4:6);
dsm_a = norm(vv_capture-V_relsoi);
    

%% Return to Earth - Definition of the launch and arrival window----------------------------- 
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

%% PorkChop analysis
[mePCdata] = porkChopDatav3(planet_2,planet_1,medep_dates,mearr_dates);
medv1_dMAX=3;
porkChopPlot(mePCdata,1,0,1,medv1_dMAX);

%% Optimization throught fmincon-------------------------------------------
lb = [medate_depi_mjd2000 medate_arri_mjd2000];
ub = [medate_depf_mjd2000 medate_arrf_mjd2000];

%x0= [sol_ga(1) sol_ga(2)];
x0 = [10496 10832;
     11261 11592;
     12101 12332];
options_fmincon = optimoptions('fmincon','UseParallel',true);

[m,n] = size(x0);
medates = zeros(m,n);
meDV = zeros(m,n);
for i = 1:size(x0,1)
medates(i,:) = fmincon(@(X)ffdv(X,2,planet_2,planet_1),x0(i,:),...
                      A,b,Aeq,beq,lb,ub,@(X)dv1conme(X),options_fmincon);
[medv_d,medv_a] = evaldv (medates(i,:),planet_2,planet_1);
meDV(i,:) = [medv_d medv_a];
end

%% Transfer orbit definition------------------------------------------------
DAYd = medates(1,1);          % Selected departure date [mjd2000]
DAYa = medates(1,2);          % Selected arrival date   [mjd2000]
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
[~,~,~,~,V1,~] = lambertMR(state_Md(1:3),state_Ea(1:3),tof_sol,mu_sun,0,0,0,0);

% Transfer orbit computation
state0 = [state_Md(1:3) V1];        % Initial state vector
kep_T = car2kep(state0,mu_sun);     % Keplerian elements of the transfer orbit
T_T = period(kep_T,mu_sun);         % Transfer orbit period

[tv_T,state_T] = ode113(@(t,y)dyn_2BP(t,y,mu_sun),[0,T_T],state0,options);


%% Escape patching
[tv_Mt,state_Mt]=ode113(@(t,y)dyn_2BP(t,y,mu_sun),0:60:tof_sol,state_Md,options);
[tv_t,state_tt]=ode113(@(t,y)dyn_2BP(t,y,mu_sun),tv_Mt,state0,options);
[rr_escape,vv_escape,date_escape] = patching(state_Mt,state_tt,tv_t,DAYd,mars_rsoi);

%% Lambert arc recalculation
kep_Msoi = uplanet(date_escape,4);
state_Msoi = kep2car(kep_Msoi, mu_sun);

[~,~,~,~,V1,~] = lambertMR(state_Msoi(1:3)+rr_escape,state_Ea(1:3),(DAYa-date_escape)*3600*24,mu_sun,0,0,0,0);
V_relsoi=V1-state_Msoi(4:6);
dsm_d = norm(vv_escape-V_relsoi)

%% Escape manoeuvres
rp_esc = kepf(1);
vv_esc_u = vv_escape/norm(vv_escape);
vesc = norm(vv_escape);
Delta = mu_mars/vesc^2*((1+vesc^2*rp_esc/mu_mars)^2-1)^0.5;

% Bplane definition for escaping
%T = 
kep_he = car2kep([rr_escape vv_escape],mu_mars);
    
    
    
    
% figure()
% hold on
% grid on
% plot(theta*180/pi,dvtot)
% dv_m = [];
% t_m = [];
% rpH = [];
% Delta = 2*mr_mars:100:3*mr_mars;
% theta = [-5:0.02:10]*pi/180;
% lb = [Delta(1) theta(1)];
% ub = [Delta(end) theta(end)];
% 
% %man_ga = ga(@(x)manopt(x,T,R,mars_rsoi,vv_capture, mu_mars, kep_p, date_capture),...
%  %   2,A,b,Aeq,beq,lb,ub,@(x)mancon(x,T,R,mars_rsoi,vv_capture, mu_mars, kep_p, date_capture),options);
% 
% sol = fmincon(@(x)manopt(x,T,R,mars_rsoi,vv_capture, mu_mars, kep_p, date_capture),[7924 theta0],...
%                       A,b,Aeq,beq,lb,ub,@(x)mancon(x,T,R,mars_rsoi,vv_capture, mu_mars, kep_p, date_capture),options_fmincon);
% 
% 
% [DVs, DVtot, T, Ttot, rp_H] = manopt(sol, T, R, rsoi, vvsoi, mu, kep_p, date_soi);
% 


% %%
% Delta = 2*mr_mars:100:3*mr_mars;
% dv_m=[];
% parfor i = 1:length(Delta)
% 
%     [dv_m(i)] = manopt([Delta(i) 23*pi/180],T,R, mars_rsoi, vv_capture, mu_mars, kep_p, date_capture)
% end
% 
% figure()
% hold on
% grid on
% plot(Delta,dv_m, 'r')

% figure()
% hold on
% grid on
% plot(Delta,t_m, 'b')
% 
% figure()
% hold on
% grid on
% plot(Delta,rpH, 'g')

% gamma = asin(norm(Delta)/mars_rsoi);
% rrpar = mars_rsoi*cos(gamma)*(-S);
% r_in = rrpar+Delta;
% state_soi_in = [r_in vv_capture];
% kep_hyp = car2kep(state_soi_in, mu_mars);
% rp_hyp = kep_hyp(1)*(1-kep_hyp(2));
% P_hyp = kep_hyp(1)*(1-kep_hyp(2)^2);
% ta_soi = acos((P_hyp/norm(r_in)-1)/kep_hyp(2));
% dt_hyp = dtheta2dt (kep_hyp,mu_mars,ta_soi);
% [~,state_hypsoi]=ode113(@(t,y)dyn_2BP(t,y,mu_mars),0:dt_hyp,state_soi_in,options);
% rrp_H = state_hypsoi(end,1:3);
% vvp_hyp = state_hypsoi(end,4:6);
% vvp_hyp_unit = vvp_hyp/norm(vvp_hyp);   % Direction of the velocity at pericenter
% 
% %% Defining the elliptical orbit
% e_el = 0.9;
% rrp_e = rrp_H;
% rp_e = rp_hyp;
% ra_e = -rp_e*(1+e_el)/(e_el-1);
% a_e = (rp_e+ra_e)/2;
% P_e = a_e*(1-e_el^2);
% vp_e = sqrt(mu_mars/P_e)*(1+e_el);
% dv_hyp_e = vp_e-norm(vvp_hyp);
% dvv_hyp_e = abs(dv_hyp_e)*(-vvp_hyp_unit);
% vvp_e = vvp_hyp + dvv_hyp_e;
% kep_e = car2kep([rrp_e vvp_e], mu_mars);
% ra_e1 = kep_e(1)*(1+kep_e(2));
% T_e = period(kep_e,mu_mars);
% [~,state_e]=ode113(@(t,y)dyn_2BP(t,y,mu_mars),0:T_e,[rrp_e vvp_e],options);
% hh = cross(rrp_e, vvp_e);
% hh = hh/norm(hh);
% kep_e(6)=pi;
% Sa1 = kep2car(kep_e,mu_mars);
% rra_e1 = Sa1(1:3);
% vva_e1 = Sa1(4:6);
% vva_e1_unit = vva_e1/norm(vva_e1);
% va_e1 = norm(vva_e1);
% 
% % Change plane
% rra_e2 = rra_e1;
% ra_e2 = ra_e1;
% rp_e2 = kep_p(1);
% e2 = (ra_e2-rp_e2)/(ra_e2+rp_e2);
% a2 = (ra_e2+rp_e2)/2;
% P_e2 = a2*(1-e2^2);
% va_e2 = sqrt(mu_mars/P_e2)*(1-e2);
% delta_i = kep_p(3)-kep_e(3);
% dv_cp = sqrt(va_e1^2+va_e2^2-2*va_e1*va_e2*cos(delta_i));
% vvh = va_e2*sin(delta_i)*-hh;
% vvpar = (va_e2*cos(delta_i)-va_e1)*vva_e1_unit;
% dvv_cp1 = vvh+vvpar;
% 
% beta = pi-asin(va_e2*sin(delta_i)/dv_cp);
% alfa = pi-beta;
% dv_par = dv_cp*cos(alfa)*vva_e1_unit;
% dv_norm = dv_cp*sin(alfa)*-hh;
% dvv_cp= dv_par+dv_norm;
% vva_e2 = vva_e1+dvv_cp;
% kep_e2 = car2kep([rra_e2 vva_e2], mu_mars);
% 
% 
% dOM = kep_p(4)-kep_e2(4);
% alf = acos(cos(kep_e2(3))^2 + sin(kep_e2(3))^2*cos(dOM));
% u1 = asin(sin(dOM)*sin(kep_e2(3))/sin(alf));
% u2 = u1;
% ta1 = u1-kep_e2(5);
% ta2 = ta1;
% w2 = u2-ta2;
% vtta = sqrt(mu_mars/P_e2)*(1+kep_e2(2)*cos(kep_e2(6)));
% dv_cOM = 2*vtta*sin(alf/2);
% 
% % Circularize
% vc = sqrt(mu_mars/rp_e2);
% vP_e2 = sqrt(mu_mars/P_e2)*(1+e2);
% dvcirc = vc-vP_e2;
