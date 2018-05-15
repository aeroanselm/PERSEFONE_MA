function [s, Ls] = mseasons(theta_degrees)

p = 251;
Ls = theta_degrees + p;

if Ls > 360
    Ls = Ls - 360;
end


if         (Ls >= 0 && Ls <= 90)
        s = 'Spring';
    elseif (Ls > 90 && Ls <= 180)
        s = 'Summer';
    elseif (Ls > 180 && Ls <= 270)
        s = 'Fall';
    elseif (Ls > 270 && Ls <= 360)
        s = 'Winter';
end
return