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
phobos.Vertices = [pvx, pvy, pvz];x0 = s0(:,1);
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
zlabel('\textbf{[km]}','interprex0 = s0(:,1);
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
zlabel('\textbf{[km]}','interprex0 = s0(:,1);
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
zlabel('\textbf{[km]}','interprex0 = s0(:,1);
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
zlabel('\textbf{[km]}','interpreter','Latex');ter','Latex');ter','Latex');ter','Latex');er','Latex');

patch('Vertices',phobos.Vertices,'Faces',phobos.Faces,'FaceColor',[0.620,0.620,0.620],'EdgeAlpha',0);
light('Position',[1000 0 10000],'Style','local');
ax = gca();
ax.DataAspectRatio = [1,1,1];

xlabel('\textbf{[km]}','interpreter','Latex');
ylabel('\textbf{[km]}','interpreter','Latex');
zlabel('\textbf{[km]}','interpreter','Latex');