function dt = dtheta2dt (kep,mu,dtheta)

% dtheta2dt.m - delta_theta to coresponding delta_time.
% 
% TYPE: 
%   Function.
%
% PROTOTYPE:
%   dt = dtheta2dt (kep,mu,dtheta)
%
% DESCRIPTION:
%   It computes the time needed to travel from the initial position and the
%   final position, which are given.   
%
%   Note:   The function is written to works in a progressive way, the
%           spacecraft is moving FORWARD.
%
%   See also:   dt2dtheta.m
%
% INPUT:
%   kep[6]      Keplerian values of the orbit.
%   mu[1]       Planetary constant of the planet (mu = mass * G) [L^3/T^2].          
%   dtheta[1]   Angle between initial and final position, MUST be >0 [rad].
%
% OUTPUT:
%   dt[1]       Time interval [s]
%
% NOTATION:
%   kep[6] is a vector made as kep=[a,e,i,Om,w,theta].
%   (if e=1 -> kep=[h,e=1,i,Om,w,theta]).
%   state[6] is a vector made as state=[x,y,z,vx,vy,vz].
%   xx is a vector.
%   x is the magnitude of vector xx.
%
% CALLED FUNCTIONS:
%   period.m
%
% FUTURE DEVELOPMENT:
%   Work in progress.
%
% --------------------------------------------------------------------------
% ***** BEGIN GPL LICENSE BLOCK *****
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
% ***** END GPL LICENCE BLOCK *****
% --------------------------------------------------------------------------
%
% AUTHORS:
%
%   Name: ALESSANDRO MARIA 
%   Surname: MASSERINI
%   ID number: 808451
%   Contact: alessandro.masserini@mail.polimi.it
%
%   Name: GIANMARIO  
%   Surname: MERISIO
%   ID number: 874400
%   Contact: gianmario.merisio@mail.polimi.it
%
%   Name: GIOVANNI ANTONIO
%   Surname: ZANOTTI
%   ID number: 876655
%   Contact: giovanni3.zanotti@mail.polimi.it
%
%   Course: Space Engineering
%   Department: DAER
%   University: Politecnico di Milano
%   Class: Orbital Mechanics 
%   Creation: 17/12/2016
%
% CHANGELOG:
%   01/01/2017, GIANMARIO MERISIO, added cases e=1 and e>1;
%
% -------------------------------------------------------------------------


% CHECK
if dtheta < 0
    error('WARNING: dtheta < 0');
end

%sinE = (sqrt(1-e^2))sentheta/1+ecostheta
%cosE = (e+costheta)/(1+ecostheta)

e = kep(2);
if e<1                  % Ellipse and Circle
    theta1 = kep(6);
    theta2=theta1+dtheta;
    n=0;
    while theta2>=2*pi
        n=n+1;
        theta2=theta2-2*pi;
    end

    senE_1 =(sqrt(1-e*e))*sin(theta1)./(1+e*cos(theta1)); 
    cosE_1 =  (e+cos(theta1))./(1+e*cos(theta1)); 
    senE_2 = (sqrt(1-e*e))*sin(theta2)./(1+e*cos(theta2));
    cosE_2 = (e+cos(theta2))./(1+e*cos(theta2));
    E_1 = atan2(senE_1,cosE_1);
    E_2 = atan2(senE_2,cosE_2);

    T = period (kep,mu);

    t_1 = T/(2*pi)*(E_1 - e*senE_1);
    t_2 = T/(2*pi)*(E_2 - e*senE_2);

    if t_1<0
        t_1=T+t_1;
    end

    if t_2<0
        t_2=T+t_2;
    end

    dt=t_2+n*T-t_1;
    
elseif e>1               % Hyperbola
    theta1 = kep(6);
    if theta1>pi
        theta1=theta1-2*pi;
    end
    theta2=theta1+dtheta;
    if theta1>theta2
        error(['DANGER (hyperbola): theta1>theta2']);
    end
    
    a = kep(1);

    H_1=2*atanh(sqrt((e-1)./(e+1)).*tan(theta1/2));
    H_2=2*atanh(sqrt((e-1)./(e+1)).*tan(theta2/2));


    t_1 = sqrt(-a^3/mu).*(e*sinh(H_1)-H_1);
    t_2 = sqrt(-a^3/mu).*(e*sinh(H_2)-H_2);

    dt = t_2 - t_1;
    
    if dt<0
        error('dt lower than 0, probably you are given wrong dtheta');
    end
    
elseif e==1              % Parabola
    theta1 = kep(6);
    if theta1>pi
        theta1=theta1-2*pi;
    end
    theta2=theta1+dtheta;
    if theta1>theta2
        error(['DANGER (parabola): theta1>theta2']);
    end
    
    h=kep(1);
    p=h^2/mu;
    
    D_1 = tan(theta1/2);
    D_2 = tan(theta2/2);
    
    t_1 = sqrt(p^3/mu).*0.5*(D_1+D_1^3/3);
    t_2 = sqrt(p^3/mu).*0.5*(D_2+D_2^3/3);
    
    dt = t_2 - t_1;
    if dt<0
        error('dt lower than 0, probably you are given wrong dtheta');
    end
    
else
    error('Eccentricity e<0');
end

return
