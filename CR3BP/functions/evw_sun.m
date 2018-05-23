function [tv,visibility_sun,visibility_mars] = evw_sun(state1, state2, tof,statep)
%This routine is used to verify the time instants in which body2 can be
%seen by body1 [Earth]. If flag == 0 only the presence of the Sun is
%considered, otherwise a secondary planet is also considered. tof can be a
%scalar or a vector indicating the time of flight that has to be
%integrated using the fcn ode113. If a secondary planet is considered the
%state2 is relative to the secondary planet (Mars).

mu_sun = astroConstants(4);
mu_mars = astroConstants(14);
AU = astroConstants(2);
h_sun = AU*sin(10*pi/180);
mr_mars = astroConstants(24);
h_mars = mr_mars + 100;

s01 = state1;
s02 = state2;

[n,m] = size(tof);
if n == 1 && m ==1
    tof = [0 tof];
end

options = odeset('Reltol',1e-13,'AbsTol',1e-14);

if nargin < 4
    [tv,s1]=ode113(@(t,y)dyn_2BP(t,y,mu_sun),tof,s01,options);
    [~,s2]=ode113(@(t,y)dyn_2BP(t,y,mu_sun),tv,s02,options);

    [n,~] = size(s1);
    visibility_sun = zeros(n,1);
    

for i = 1:n
    AA = s1(i,1:3);
    BB = s2(i,1:3);
    CC = BB - AA;
    A = norm(AA);
    B = norm(BB);
    C = norm(CC);
    p = (A + B + C)/2;
    area = sqrt(p*(p - A)*(p - B)*(p - C));
   
    H = 2*area/C;
    if H <= h_sun
        visibility_sun(i) = 0;
    else
        visibility_sun(i) = 1;
    end
end
else
    s0p = statep;
    [tv,s1]=ode113(@(t,y)dyn_2BP(t,y,mu_sun),tof,s01,options);
    [~,sp]=ode113(@(t,y)dyn_2BP(t,y,mu_sun),tv,s0p,options);
    [~,s2]=ode113(@(t,y)dyn_2BP(t,y,mu_mars),tv,s02,options);
    [n,~] = size(s1);
    visibility_sun = zeros(n,1);
    visibility_mars = zeros(n,1);
    
    for i = 1:n
    AA = sp(i,1:3) - s1(i,1:3);
    BB = s2(i,1:3);
    CC = AA + BB;
    A = norm(AA);
    B = norm(BB);
    C = norm(CC);
    p = (A + B + C)/2;
    area = sqrt(p*(p - A)*(p - B)*(p - C));
    H = 2*area/C;
    if H <= h_mars && C > A
        visibility_mars(i) = 0;
    else
        visibility_mars(i) = 1;
    end
    end
    
    for i = 1:n
    AA = s1(i,1:3);
    BB = sp(i,1:3) + s2(i,1:3);
    CC = BB - AA;
    A = norm(AA);
    B = norm(BB);
    C = norm(CC);
    p = (A + B + C)/2;
    area = sqrt(p*(p - A)*(p - B)*(p - C));
    H = 2*area/C;
    if H <= h_sun
        visibility_sun(i) = 0;
    else
        visibility_sun(i) = 1;
    end
    end
    
end

    

return
    
