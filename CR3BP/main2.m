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
d = 9375;
M = m_mars + m_phobos;
p1 = m_mars/M;
p2 = m_phobos/M;

x1 = -p2*d;
x2 = p1*d;


y0 = [9360 3 0 0 20 0 3600/2];
x = fsolve(@(y)correction_dim(m_mars,m_phobos,mu_mars,mu_phobos,d,y),y0);
[t,state]=ode113(@(t,x)dynCR3BPdim(t,x,d,m_mars,m_phobos,mu_mars,mu_phobos),[0 80*x(7)],x,options);
dv = state(end,4:6) - state(1,4:6);
dy = state(end,1:3) - state(1,1:3);

% Libration points
[L,~] = libpoints(m_mars,m_phobos,d);

figure()
hold on
grid on
axis equal
drawPlanet('Mars',[x1,0,0],gca,1);
drawPlanet('Phobos',[x2,0,0],gca,1);
plot3(state(:,1),state(:,2),state(:,3))
