clear
close all
clc

addpath(genpath('./functions'));
options = odeset('Reltol',1e-13,'AbsTol',1e-16);

G = astroConstants(1);                  
mu_mars = astroConstants(14);
m_mars = mu_mars/G;
mu_phobos = 0.0007112;
m_phobos = mu_phobos/G;

mu = m_phobos/(m_mars + m_phobos);

y0 = [4.004952718E+03;...
      8.517437775E+03;...
     -5.589719112E+01;...
     -1.908368039E+00;...
      9.321397618E-01;...
     -3.865928247E-02;...
     3600*0.5];

foptions = optimoptions('fsolve','FunctionTolerance',1e-13,'OptimalityTolerance',1e-14,'StepTolerance',1e-14);
x = fsolve(@(y)correction(mu,y),y0,foptions);
[t,state]=ode113(@(t,x)dyn_CR3BP1(t,x,mu),[0 2*x(7)],x,options);
dv = state(end,4:6) - state(1,4:6);
dy = state(end,1:3) - state(1,1:3);

figure()
hold on
grid on
axis equal
plot3(state(:,1),state(:,2),state(:,3))