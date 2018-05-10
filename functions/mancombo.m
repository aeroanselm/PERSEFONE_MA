function [dvv, dv, dt,stateall, rp_h, kepf,sol] = mancombo(state,mu,i_f,OM_f,e_heo,rp_f)

options = odeset('Reltol',1e-13,'AbsTol',1e-14);

% option 1
state_h = state;
rp2 = rp_f;
kep_h = car2kep(state_h,mu);
if kep_h<1
    error('The state at the entry of the soi is not an hyperbola')
end

[dvv_h_heo, dv_h_heo, dt_h_heo, kep_heo, rp_h] = hyp2heo(state_h,mu,e_heo);
[dvv_ps, dv_ps, dt_ps, kep_heo2] = changeps_apo(kep_heo,rp2,mu,i_f);
state02 = kep2car(kep_heo,mu);

dOM = OM_f - kep_heo2(4);
[dvv_OM, dv_OM, dt_OM, kep2] = changeOM(kep_heo2, mu, dOM);
state03 = kep2car(kep_heo2,mu);

[dvv_c, dv_c, dt_c, kepf] = orbit_circ(kep2, mu);
state04 = kep2car(kep2,mu);
T = period(kepf,mu);
state05 = kep2car(kepf,mu);

[~,state_1]=ode113(@(t,y)dyn_2BP(t,y,mu),0:1:dt_h_heo,state_h,options);
[~,state_2]=ode113(@(t,y)dyn_2BP(t,y,mu),0:1:dt_ps,state02,options);
[~,state_3]=ode113(@(t,y)dyn_2BP(t,y,mu),0:1:dt_OM,state03,options);
[~,state_4]=ode113(@(t,y)dyn_2BP(t,y,mu),0:1:dt_c,state04,options);
[~,state_5]=ode113(@(t,y)dyn_2BP(t,y,mu),0:1:T,state05,options);

stateall1 = [state_1;state_2;state_3;state_4;state_5];

dvv1 = [dvv_h_heo; dvv_ps; dvv_OM; dvv_c];
dv1 = [dv_h_heo; dv_ps; dv_OM; dv_c];
dt1 = [dt_h_heo; dt_ps; dt_OM; dt_c];

    
% Option 2
[dvv_h_heo, dv_h_heo, dt_h_heo, kep_heo] = hyp2heo(state_h,mu,e_heo);
[dvv_ps, dv_ps, dt_ps, kep_heo2] = changeps_apo(kep_heo,rp2,mu,i_f);
state02 = kep2car(kep_heo,mu);

[dvv_c, dv_c, dt_c, kep_c] = orbit_circ(kep_heo2, mu);
state03 = kep2car(kep_heo2,mu);

dOM = OM_f - kep_c(4);
[dvv_OM, dv_OM, dt_OM, kepf] = changeOM(kep_c, mu, dOM);
state04 = kep2car(kep_c,mu);
state05 = kep2car(kepf,mu);
T = period(kepf,mu);
dvv2 = [dvv_h_heo; dvv_ps; dvv_OM; dvv_c];
dv2 = [dv_h_heo; dv_ps; dv_OM; dv_c];
dt2 = [dt_h_heo; dt_ps; dt_OM; dt_c];

[~,state_12]=ode113(@(t,y)dyn_2BP(t,y,mu),0:1:dt_h_heo,state_h,options);
[~,state_22]=ode113(@(t,y)dyn_2BP(t,y,mu),0:1:dt_ps,state02,options);
[~,state_32]=ode113(@(t,y)dyn_2BP(t,y,mu),0:1:dt_c,state03,options);
[~,state_42]=ode113(@(t,y)dyn_2BP(t,y,mu),0:1:dt_OM,state04,options);
[~,state_52]=ode113(@(t,y)dyn_2BP(t,y,mu),0:1:T,state05,options);
stateall2 = [state_12;state_22;state_32;state_42;state_52];

if sum(dv1) < sum(dv2)
    dvv = dvv1;
    dv = dv1;
    dt = dt1;
    stateall = stateall1;
    sol = 1;
else
    dvv = dvv2;
    dv = dv2;
    dt = dt2;
    stateall = stateall2;
    sol = 2;
end


