clear
close all
clc
addpath(genpath('./functions'));
options = odeset('Reltol',1e-13,'AbsTol',1e-16);
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

x1 = (m_phobos/m_mars)/(1 + m_phobos/m_mars);
x2 = m_mars/m_phobos*x1;

mu = m_phobos/(m_mars + m_phobos);
[n,m] = size(x0);

figure()
hold on;
grid on;
axis equal;
title('\textbf{DROs around Phobos}','interpreter','Latex');
%drawPlanet('Phobos',[x2*DU 0 0],gca,1);

xlabel('\textbf{[km]}','interpreter','Latex');
ylabel('\textbf{[km]}','interpreter','Latex');
for i = 15
    s0 = [(x0(i)+x2) 0 0 0 Vy0(i) 0];
    [~,state]=ode113(@(t,x)dynCR3BPpert(t,x,mu),[0 8*pi],s0,options);
    xs = state(:,1:3)*DU;
    vs = state(:,4:6)*DU/TU;
    state = [xs vs];
    phobos = stlread('Phobos_200k.stl');
    pvx = phobos.Vertices(:,1)+x2*DU;
    pvy = phobos.Vertices(:,2);
    pvz = phobos.Vertices(:,3);
    phobos.Vertices = [pvx, pvy, pvz];
    patch('Vertices',phobos.Vertices,'Faces',phobos.Faces,'FaceColor',[0.455,0.380,0.365],'EdgeAlpha',0);
    light 
    ax = gca();
    ax.DataAspectRatio = [1,1,1];

    plot3(state(:,1),state(:,2),state(:,3));
end




