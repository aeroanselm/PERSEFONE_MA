
clear
close all
clc

%% Option #1

A0 = 15;
B0 = -15:15;
n = sqrt(42828/9379^3);
t = 0:2*pi/n;
alfa = 0;
beta = 0;
x_off = 0;
y_off = 0;

figure()
hold on
[X,Y,Z] = sphere(1e2);
X = 11*X; Y = 11*Y; Z = 11*Z;
fvc = surf2patch(X,Y,Z);
patch('Vertices',fvc.vertices,'Faces',fvc.faces);


for i = 1:length(B0)
x = A0*cos(n*t + alfa) + x_off;
y = -2*A0*sin(n*t + alfa) -3/2*n*t*x_off + y_off;
z = B0(i)*cos(n*t + beta);

% xd = -A0*sin(n*t+alfa);
% yd = -2*n*A0*cos(n*t+alfa);
% zd = -n*B0(i)*sin(n*t+beta);

% S(i,:) = 
plot3(x,y,z)
end
ax = gca();
ax.DataAspectRatio = [1,1,1];

vv = [];
for k = 1:length(B0)
t = pi/(2*n);
xd = -A0*sin(n*t+alfa);
yd = -2*n*A0*cos(n*t+alfa);
zd = -n*B0(k)*sin(n*t+beta);
vv(k,:) = [xd yd zd];
end
dv_tot = norm(vv(2,:)-vv(1,:))*length(B0);

%% Option #2

A0 = 15;
B0 = -15:15;
n = sqrt(42828/9379^3);
t = 0:2*pi/n;
alfa = 0;
beta = 0;
x_off = 0;
y_off = 0;
