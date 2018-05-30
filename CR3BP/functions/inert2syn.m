function [synodic] = inert2syn(state,t)
% This routine makes the transformation from inertial to synodic (rotating)
% coordinates. OM is set to 1; state is a 3x1 vector.

[n,m] = size(state);
if n~=3 && m~=1
    error('state must be a 3x1 vector');
end

R = [cos(t) sin(t) 0;...
    -sin(t) cos(t) 0;...
        0     0    1];
    
synodic = R*state;
return