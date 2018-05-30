function [SR] = rotmodel(state,angle)

R = [cos(angle) sin(angle) 0;...
    -sin(angle) cos(angle) 0;...
        0           0        1];
[n,m] = size(state);
SR = zeros(n,m);
for i = 1:n
    S = state(i,:);
    S = S';
    SR_temp = R*S;
    SR_temp = SR_temp';
    
    SR(i,:) = SR_temp;
end
return