function [visibility_dates] = evw (mr_body1,state_body1, state_body2, dates,flag_plot,altitude)

if length(dates)~=size(state_body1,1) || ... 
   size(state_body2,1)~=size(state_body1,1) || ...
   size(state_body2,1)~=length(dates)
disp("ERROR: matrix dimesions don't match")
return
end

if nargin<6
    altitude = 100;
end
    

rr1 = [state_body1(:,1), state_body1(:,2), state_body1(:,3)];
rr2 = [state_body2(:,1), state_body2(:,2), state_body2(:,3)];
rr12 = rr2-rr1;

visibility_dates = [];
for k = 1:length(dates)
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
        visibility_dates(k)=dates(k);
    end
end

if flag_plot==1
    bool = [];
    for m = 1:length(visibility_dates)
        if visibility_dates(m)~=NaN
            bool(m)=1;
        else
            bool(m)=0;
        end
    end
    figure()
    hold on
    grid on
    plot(visibility_dates,bool)
end
    
