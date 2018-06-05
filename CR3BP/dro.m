clear
close all
clc
addpath(genpath('./functions'));
options = odeset('Reltol',1e-13,'AbsTol',1e-14);
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
title('\textbf{DROs around Phobos}','interpreter','Latex');
phobos = stlread('Phobos_200k.stl');
angle = 1/2*pi;
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
for i = 20
    s0 = [(x0(i)+x2) 0 0 0 Vy0(i) 0];
    t = 13*pi;
    [tv,state]=ode113(@(t,x)dynCR3BPpert(t,x,mu),[0:0.001:1*pi],s0,options);
    tv = tv*TU;
%     save('tv.mat','tv');
    xs = state(:,1:3)*DU;
    vs = state(:,4:6)*DU/TU;
    state_dro = [xs vs];
%     save('state.mat','state');
    T = t*TU/3600;
    plot3(state_dro(:,1),state_dro(:,2),state_dro(:,3),'Linewidth',1.5);
end

%% RV
d = (x2 - x1);
xsc = 9347/DU;
vsc = sqrt(mu_mars/xsc/DU);
% s0 = [xsc 0 0 0 vsc*TU/DU 0];


OM = sqrt(MU/DU^3);
vsyn = OM*xsc*DU;
dv_a = (vsc - vsyn)*TU/DU;
p = 2*pi*sqrt(9354^3/mu_mars);
p_a = p/TU;

s0 = [xsc*cos(1*pi/180) -xsc*sin(1*pi/180)  0 dv_a*sin(1*pi/180)  dv_a*cos(1*pi/180)  0];
%s0 = [9358.06141509580/DU -13.0485979370938/DU 0 0.00368977660351285/DU*TU 0.00887246736567869/DU*TU 0];
w = dv_a / xsc;
period = 2*pi/w;
%%
%p_sc = [p_sc(1:3)/DU p_sc(4:6)/DU*TU];
[~,state]=ode113(@(t,x)dynCR3BPpert(t,x,mu),[0:0.001:2*pi],s0,options);
xs = state(:,1:3)*DU;
vs = state(:,4:6)*DU/TU;
state_sc = [xs vs];
T = t*TU/3600;


[dv_dro,p_sc,i] = drocap(state_dro, state_sc);
plot3(state_sc(1:i,1),state_sc(1:i,2),state_sc(1:i,3));%,'Linewidth',1.5);

plot3(p_sc(1),p_sc(2),p_sc(3),'o');


p_sc = [p_sc(1:3)/DU (p_sc(4:6) - dv_dro)/DU*TU];
[~,state]=ode113(@(t,x)dynCR3BPpert(t,x,mu),[0:0.1:130*pi],p_sc,options);
xs = state(:,1:3)*DU;
vs = state(:,4:6)*DU/TU;
state_sc = [xs vs];
T = t*TU/3600;
comet3(state_sc(:,1),state_sc(:,2),state_sc(:,3));%,'y','Linewidth',1.5);
 
function [dv_dro,p_sc,i] = drocap(state_dro, state_sc)
    tol = 0.02;
    [n,m] = size(state_sc);
    
    for i = 1:n
        rdro = [state_dro(:,1) state_dro(:,2) state_dro(:,3)];
        rsc = [state_sc(i,1) state_sc(i,2) state_sc(i,3)];
        d = rsc - rdro;
        for k = 1:size(state_dro,1)
            error = norm(d(k,:));
            if error <= tol
                p_dro = [state_dro(k,1) state_dro(k,2) state_dro(k,3)];
                p_sc = [state_sc(i,1) state_sc(i,2) state_sc(i,3) state_sc(i,4) state_sc(i,5) state_sc(i,6)];
                vv_dro = [state_dro(k,4) state_dro(k,5) state_dro(k,6)];
                vv_sc = [state_sc(i,4) state_sc(i,5) state_sc(i,6)];
                dv_dro = vv_sc - vv_dro;
                return
            end
        
        end
    end
end



