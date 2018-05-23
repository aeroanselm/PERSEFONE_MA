function [visibility_dates] = evw (mr_body1,state_body1, state_body2, dates,flag_plot,altitude)


if length(dates)~=size(state_body1,1) || ... 
   size(state_body2,1)~=size(state_body1,1) || ...
   size(state_body2,1)~=length(dates)
error("Matrix dimesions don't match")
return
end

if nargin<6
    altitude = 100;
end
    

rr1 = [state_body1(:,1), state_body1(:,2), state_body1(:,3)];
rr2 = [state_body2(:,1), state_body2(:,2), state_body2(:,3)];
rr12 = rr2-rr1;

visibility_dates = [];
mu_sun = astroConstants(4);             % Sun planetary constant
mr_sun = astroConstants(3);

for k = 1:length(dates)
    kep_E = uplanet(dates(k),3);
    kep_M = uplanet(dates(k),4);
    SE = kep2car(kep_E,mu_sun);
    SM = kep2car(kep_M,mu_sun);
    RE = SE(1:3);   re = norm(RE);
    RM = SM(1:3);   rm = norm(RM);
    REM = RM-RE;    rem = norm(REM);
    ps = (re + rm + rem)/2;
    As = sqrt(ps*(ps-re)*(ps-rm)*(ps-rem));
    hs = 2*As/rem;
    
    if hs <= (mr_sun+100000)
        visibility_dates(k)=NaN;
    else
    
        r1 = [rr1(k,1) rr1(k,2) rr1(k,3)];
        r2 = [rr2(k,1) rr2(k,2) rr2(k,3)];
        r12 = [rr12(k,1) rr12(k,2) rr12(k,3)];
        AC = norm(r1);
        AB = norm(r2);
        CB = norm(r12);
        p = (AC+AB+CB)/2;
        A = sqrt(p*(p-AC)*(p-AB)*(p-CB));
        h = 2*A/AB;
    
        if h<=(mr_body1+altitude)
            visibility_dates(k)=NaN;
        else
            visibility_dates(k)= dates(k);
        end
    end
end

if flag_plot==1
    bool = [];
    for m = 1:length(visibility_dates)
        if isnan(visibility_dates(m))
            bool(m)=0;
        else
            bool(m)=1;
        end
    end
    figure()
    hold on
    grid on
    date1 = mjd20002date(visibility_dates(1));
    date2 = mjd20002date(visibility_dates(end));
    plot(visibility_dates,bool,'linewidth',2)
    title('Earth visibility window from 04/10/2025 21:35 to 06/10/2025 21:35')
%    subtitle ('Eclipses duration ~35 min Exposure duration ~3h')
end
    
