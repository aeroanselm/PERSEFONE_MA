clear
close all
clc
%%
addpath(genpath('./functions'));
options = odeset('Reltol',1e-13,'AbsTol',1e-14);
foptions = optimoptions('fsolve','FunctionTolerance',1e-13,'OptimalityTolerance',1e-14,'StepTolerance',1e-14);
load('IC.mat');

x0 = s0(:,1);
Vy0 = s0(:,2);
DU = 9379.11746340512;
TU = 4389.11709709003;

G = astroConstants(1);                  
mu_mars = astroConstants(14);
m_mars = mu_mars/G;
mu_phobos = 0.0007112;
m_phobos = mu_phobos/G;
p1 = m_mars/(m_mars + m_phobos);
p2 = m_phobos/(m_mars + m_phobos);
M = m_mars + m_phobos;
MU = G*M;
x1 = -p2;
x2 = p1;

mu = m_phobos/(m_mars + m_phobos);
[n,m] = size(x0);

figure()
hold on;
grid on;
axis equal;
title('\textbf{DRO at 10 km from Phobos surface}','interpreter','Latex');
phobos = stlread('Phobos_200k.stl');
angle = pi;
phobos.Vertices = rotmodel(phobos.Vertices,angle);
pvx = phobos.Vertices(:,1)+x2*DU;
pvy = phobos.Vertices(:,2);
pvz = phobos.Vertices(:,3);
phobos.Vertices = [pvx, pvy, pvz];

patch('Vertices',phobos.Vertices,'Faces',phobos.Faces,'FaceColor',[0.620,0.620,0.620],'EdgeAlpha',0);
light('Position',[1000 0 10000],'Style','local');
ax = gca();
ax.DataAspectRatio = [1,1,1];

xlabel('\textbf{[km]}','interpreter','Latex');
ylabel('\textbf{[km]}','interpreter','Latex');
zlabel('\textbf{[km]}','interpreter','Latex');

i = 20;
y0 = [(x0(i)+x2) 0 0 0 Vy0(i) 0 2*0.7*pi];
%y0 = [(40/DU + x2) 0 0 0 sqrt(mu_phobos/40)/DU*TU 0 (2*pi*sqrt(40^3/mu_phobos)/TU

x = fmincon(@(y)correction_dro(mu,y,1),y0,[],[],[],[]);

%x = fsolve(@(y)correction_dro(mu,y,0),y0,foptions);
% figure()
% hold on
% grid on

[tv,state]=ode113(@(t,x)dynCR3BPpert(t,x,mu,0),[0 1*x(7)],x,options);


xs = state(:,1:3)*DU;
vs = state(:,4:6)*DU/TU;
state_dro = [xs vs];
plot3(state_dro(:,1),state_dro(:,2),state_dro(:,3));%,'Linewidth',1.5);

