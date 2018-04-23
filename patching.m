function [rr,vv,date] = patching(state1,state2,T,datei,r_soi)

rr_planet = [state1(:,1),state1(:,2),state1(:,3)];
rr_body = [state2(:,1),state2(:,2),state2(:,3)];
vv_planet = [state1(:,4),state1(:,5),state1(:,6)];
vv_body = [state2(:,4),state2(:,5),state2(:,6)];
Tdays = T/3600/24;
rr_pb = rr_body-rr_planet;
vv_pb = vv_body-vv_planet;

tol = 100;
r = [];
for k = 1:length(T)
    r_pb = norm([rr_pb(k,1), rr_pb(k,2), rr_pb(k,3)]);
    r =[r;r_pb];
    if abs(r_pb-r_soi)<tol
        rr = [rr_pb(k,1), rr_pb(k,2), rr_pb(k,3)];
        vv = [vv_pb(k,1), vv_pb(k,2), vv_pb(k,3)];
        date = datei+Tdays(k);
        return
    else
        rr = NaN;
        vv = NaN;
        date = NaN;
    end
end