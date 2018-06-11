%% Continuation algorithm for DRO in Mars-Phobos CRTBP system
clear
close all
clc
%% Data section
addpath(genpath('./functions'));
options = odeset('Reltol',1e-13,'AbsTol',1e-14);

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
T = 2*pi*sqrt(DU^3/MU);
OM = 2*pi/T;

load('IC_DRO.mat');
%% Construction of the initial guess
r2BP = 13;
x_a = (x2 - r2BP/DU);
v2BP = sqrt(mu_phobos/r2BP);
vsyn = OM*(x_a);
Dv = v2BP - vsyn;
p = 2*pi*sqrt(r2BP^3/mu_phobos);
%% Correction algorithm applied at first guess
% x0 = [x_a 0 0 0 Dv 0 p/TU];
x0 = IC(end,:);
y0 = [x0(5) x0(7)];
x = fmincon(@(y)correction_dro1(mu,x0(1),y,1),y0,[],[],[],[]);
x = [x0(1) 0 0 0 x(1) 0 x(2)];
[tv,state]=ode113(@(t,x)dynCR3BPpert(t,x,mu,0),[0 1*x(7)],x,options);
xs = state(:,1:3)*DU;
vs = state(:,4:6)*DU/TU;
state_dro = [xs vs];

%% Continuation of the family
x_step = 1/DU;
n = 2000;
IC = zeros(n,7);

for i = 1:n
    x0 = [(x(1) - x_step) 0 0 0 x(5) 0 x(7)];
    y0 = [x0(5) x0(7)];
    x = fmincon(@(y)correction_dro1(mu,x0(1),y,1),y0);
    x = [x0(1) 0 0 0 x(1) 0 x(2)];
    IC(i,:) = x;
end
%% Plot 
figure()
hold on;
grid on;
axis equal;
title('\textbf{DRO 1 to 1000 km from Phobos surface (100 km step)}','interpreter','Latex');
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
drawPlanet('Mars',[0 0 0],gca,1);
xlabel('\textbf{[km]}','interpreter','Latex');
ylabel('\textbf{[km]}','interpreter','Latex');
zlabel('\textbf{[km]}','interpreter','Latex');
%%
for i = 1:100:size(IC,1)
    [tv,state]=ode113(@(t,x)dynCR3BPpert(t,x,mu,0),[0 100*IC(i,7)],IC(i,:),options);
    xs = state(:,1:3)*DU;
    vs = state(:,4:6)*DU/TU;
    state_dro = [xs vs];
    plot3(state_dro(:,1),state_dro(:,2),state_dro(:,3));
end
