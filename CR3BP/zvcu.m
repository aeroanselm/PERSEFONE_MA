function [out] = zvcu(x,y,m1,m2,r12)
% x = x(1,:);
% y = y(:,1);

G = 6.67e-20;
p1 = m1/(m1+m2);
p2 = m2/(m1+m2);
mu = G*(m1+m2);
Om = sqrt(mu/r12^3);
mu1 = G*m1;
mu2 = G*m2;

% r1 = @(x,y) sqrt((x+p2*r12)^2+y^2);
% r2 = @(x,y) sqrt((x-p1*r12)^2+y^2);
% f = @(x,y) Om^2*(x^2+y^2)+2*mu1/r1+2*mu2/r2;
% 
% 
% C = out/(-2);

out = [];
for k = 1:length(x)
    for l = 1:length(y)
        r1 = sqrt((x(k,l)+p2*r12)^2+y(k,l)^2);
        r2 = sqrt((x(k,l)-p1*r12)^2+y(k,l)^2);
        out(k,l) = -(Om^2*(x(k,l)^2+y(k,l)^2)+2*mu1/r1+2*mu2/r2)/2;%+2*C;
        
    end
end


