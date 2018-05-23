clear all
close all
clc

%% Constants definition Earth-Moon
r12 = 384400;
m1 = 5.972e24;
m2 = 7.342e22;
M = m1+m2;
G = 6.67e-20;
p1 = m1/(m1+m2);
p2 = m2/(m1+m2);
mu = G*M;
Om = sqrt(mu/(r12^3));

%% Mars-Phobos
r12 = 9375;
m1 = 6.4185e23;
m2 = 1.07e16;
M = m1+m2;
G = 6.67e-20;
p1 = m1/(m1+m2);
p2 = m2/(m1+m2);
mu = G*M;
Om = sqrt(mu/(r12^3));

% % Synodic reference
% d1 = -p2*r12;
% d2 = p1*r12;
% 
% % Lagrange Points 
% 
% x4 = r12/2-p2*r12;
% x5 = x4;
% y4 = +sqrt(3)/2*r12;
% y5 = -sqrt(3)/2*r12;
% 
% fzeta = @(zeta) ((1-p2)./abs(zeta+p2).^3) .* (zeta+p2) + p2./abs(zeta+p2-1).^3 .* (zeta+p2-1)-zeta;
% x1 = fzero(fzeta,0.8)*r12;
% x2 = fzero(fzeta,-1)*r12;
% x3 = fzero(fzeta,1.1)*r12;
% 
% % figure()
% % hold on
% % grid on
% % plot([-2:0.01:2],fzeta([-2:0.01:2]))
% 
% ratio = 10;
% 
% 
% % Zero velocity curves
% mu1 = G*m1;
% mu2 = G*m2;
% % x = [-x3:1000:x3];
% %y = [-x3:1000:x3];
% % r1 = @(x,y) sqrt((x+p2*r12).^2+y.^2);
% % r2 = @(x,y) sqrt((x-p1*r12).^2+y.^2);
% 
% [x,y] = meshgrid([-6e5:1000:6e5],[-6e5:1000:6e5]);
% % r1 =  sqrt((x+p2*r12).^2+y.^2);
% % r2 =  sqrt((x-p1*r12).^2+y.^2);
% %%
% C = -1.6735;
% 
% %zvcu = @(x,y) Om^2.*(x.^2+y.^2)+2*mu1./r1(x,y)+2*mu2./r2(x,y) + 2*C;
% %zvcu =  Om^2.*(x.^2+y.^2)+2*mu1./r1+2*mu2./r2 + 2*C;
% %pool = parpool
% [out] = zvcu(x,y,m1,m2,r12)
% 
% 
% % figure()
% % hold on
% % grid on
% % axis equal
% % plot3(d1,0,0,'o');
% % plot3(d2,0,0,'o');
% % plot3(x1,0,0,'x');
% % plot3(x2,0,0,'x');
% % plot3(x3,0,0,'x');
% % plot3(x4,y4,0,'x');
% % plot3(x5,y5,0,'x');
% %drawPlanet('Earth',[d1 0 0],gca,ratio);
% %drawPlanet('Moon',[d2 0 0],gca,ratio);
% 
% 
% figure()
% hold on
% grid on
% axis equal
% contourf(x,y,out,[C,C])
% % plot3(d1,0,0,'o');
% % plot3(d2,0,0,'o');
% plot3(x1,0,0,'x');
% plot3(x2,0,0,'x');
% plot3(x3,0,0,'x');
% plot3(x4,y4,0,'x');
% plot3(x5,y5,0,'x');
% drawPlanet('Earth',[d1 0 0],gca,ratio);
% drawPlanet('Moon',[d2 0 0],gca,ratio);
%%
body_1 = "Earth";
body_2 = "Moon";
x = [-6e5:1000:6e5];
y =x;
C = -1.6705;
zvcv1(x,y,body_1,body_2,m1,m2,r12,C);
%%
fzeta = @(zeta) ((1-p2)./abs(zeta+p2).^3) .* (zeta+p2) + p2./abs(zeta+p2-1).^3 .* (zeta+p2-1)-zeta;
x1 = fzero(fzeta,0.8)*r12;
body_1 = 'Mars';
body_2 = 'Phobos';
x = [-1e4:10:1e4];
y =x;
C = -6.849899;
[out,d1,d2]=zvcv1(x,y,body_1,body_2,m1,m2,r12,C);

%%
[L,d]=libpoints(m1,m2,r12);

figure()
hold on
grid on
axis equal
plot(L(:,1),L(:,2),'o',d(1),0,'x',d(2),0,'x')

