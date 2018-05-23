clear
close all
clc

addpath(genpath('./functions'));
options = odeset('Reltol',1e-13,'AbsTol',1e-16);

G = astroConstants(1);                  
mu_earth = astroConstants(13);
m_earth = mu_earth/G;
mu_moon = astroConstants(20);
m_moon = mu_moon/G;
d = 384400;
M = m_earth + m_moon;
p1 = m_earth/M;
p2 = m_moon/M;

x1 = -p2*d;
x2 = p1*d;


y0 = [4e5 5 0 0 5 0 3600*8];
x = fsolve(@(y)correction_dim(m_earth,m_moon,mu_earth,mu_moon,d,y),y0);
[t,state]=ode113(@(t,x)dynCR3BPdim(t,x,d,m_earth,m_moon,mu_earth,mu_moon),[0 80*x(7)],x,options);
dv = state(end,4:6) - state(1,4:6);
dy = state(end,1:3) - state(1,1:3);

% Libration points
[L,~] = libpoints(m_earth,m_moon,d);

figure()
hold on
grid on
axis equal
drawPlanet('Earth',[x1,0,0],gca,1);
drawPlanet('Moon',[x2,0,0],gca,1);
plot3(state(:,1),state(:,2),state(:,3))
