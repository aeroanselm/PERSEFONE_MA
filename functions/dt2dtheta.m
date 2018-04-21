function dtheta = dt2dtheta (kep,mu,dt)

% dt2dtheta.m - delta_time to coresponding delta_theta.
% 
% TYPE: 
%   Function.
%
% PROTOTYPE:
%   dtheta = dt2dtheta (kep,mu,dt)
%
% DESCRIPTION:
%   It computes the angle that passes from the starting position after an
%   interval of time that is given.    
%
%   Note:   The function is written to works in a progressive way, the
%           spacecraft is moving FORWARD.
%
%   See also:   dtheta2dt.m
%
% INPUT:
%   kep[6]      Keplerian values of the orbit.
%   mu[1]       Planetary constant of the planet (mu = mass * G) [L^3/T^2]          
%   dt[1]       Time interval [s]
%
% OUTPUT:
%	dtheta[1]   Angle between initial and final position [rad]
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
%   dd/mm/yyyy, NAME SURNAME, informations
%
% -------------------------------------------------------------------------

% CHECK
if dt < 0
    error('WARNING: dt < 0');
end

%senE = (sqrt(1-e^2))sentheta/1+ecostheta
%cosE = (e+costheta)/(1+ecostheta)
e = kep(2);

if e<1                  % Ellipse and Circle
    theta1 = kep(6);
    senE_1 =(sqrt(1-e*e))*sin(theta1)./(1+e*cos(theta1)); 
    cosE_1 =  (e+cos(theta1))./(1+e*cos(theta1)); 
    E_1 = atan2(senE_1,cosE_1);

    T = period(kep,mu);

    t_1 = T/(2*pi)*(E_1 - e*senE_1);
    t_2 = t_1+dt;

    fun=@(x) t_2-T/(2*pi)*(x - e*sin(x));
    E_2=fzero(fun,0);

    theta2=2*atan(sqrt((1+e)/(1-e))*tan(E_2/2));

    while (t_2>T)
        theta2=2*pi+theta2;
        t_2=t_2-T;
    end

    dtheta=theta2-theta1;
    
elseif e>1              % Hyperbola
    theta1 = kep(6);
    if theta1>pi
        theta1=theta1-2*pi;
    end
    
    a=kep(1);
    
    H_1=2*atanh(sqrt((e-1)./(e+1)).*tan(theta1/2));
    t_1 = sqrt(-a^3/mu).*(e*sinh(H_1)-H_1);
    t_2=t_1+dt;
    
    fun=@(x) t_2-sqrt(-a^3/mu).*(e*sinh(x)-x);
    H_2=fzero(fun,0);
    
    theta2=2*atan(sqrt((e+1)/(e-1))*tanh(H_2/2));
    
    dtheta=theta2-theta1;
    
elseif e==1              % Parabola
    theta1 = kep(6);
    if theta1>pi
        theta1=theta1-2*pi;
    end
    
    h=kep(1);
    p=h^2/mu;
    
    D_1 = tan(theta_1/2);
    t_1 = sqrt(p^3/mu).*0.5*(D_1+D_1^3/3);
    t_2 = t_1+dt;
    
    fun=@(x) t_2-sqrt(p^3/mu).*0.5*(x+x.^3/3);
    D_2=fzero(fun,0);
    
    theta2=2*atan(D_2);
    
    dtheta=theta2-theta1;
    
else
        error('Eccentricity e<0');
end

  
return
