function [out,d1,d2]=zvcv1(x,y,body_1,body_2,m1,m2,r12,C)

G = 6.67e-20;
p1 = m1/(m1+m2);
p2 = m2/(m1+m2);
mu = G*(m1+m2);
Om = sqrt(mu/r12^3);
mu1 = G*m1;
mu2 = G*m2;

% Synodic reference
d1 = -p2*r12;
d2 = p1*r12;

[x,y] = meshgrid(x,y);
out = [];
for k = 1:length(x)
    for l = 1:length(y)
        r1 = sqrt((x(k,l)+p2*r12)^2+y(k,l)^2);
        r2 = sqrt((x(k,l)-p1*r12)^2+y(k,l)^2);
        out(k,l) = -(Om^2*(x(k,l)^2+y(k,l)^2)+2*mu1/r1+2*mu2/r2)/2; %=C;
        
    end
end

ratio = 1;

figure()
hold on
grid on
axis equal
contourf(x,y,out,[C,C])
drawPlanet(body_1,[d1 0 0],gca,ratio);
drawPlanet(body_2,[d2 0 0],gca,ratio);
