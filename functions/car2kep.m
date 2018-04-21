function [kep] = car2kep (state, mu)

% car2kep.m - From cartesian to keplerian coordinates.
% 
% TYPE: 
%   Function.
%
% PROTOTYPE:
%   [kep] = car2kep (state, mu)
%
% DESCRIPTION:
%   The function uses the standard algorithm to pass from cartesian 
%   coordinates to keplerian coordinates or orabital parameters.
%   
%   See also:   kep2car.m
%
% INPUT:
%	state[6]    It is the state of the spacecraft in the cartesian
%               reference frame
%   mu[1]       Planetary constant of the planet (mu = mass * G) [L^3/T^2]          
%
% OUTPUT:
%   kep[6]      Keplerian values of the orbit.
%
% NOTATION:
%   kep[6] is a vector made as kep=[a,e,i,Om,w,theta].
%   (if e=1 -> kep=[h,e=1,i,Om,w,theta]).
%   state[6] is a vector made as state=[x,y,z,vx,vy,vz].
%   xx is a vector.
%   x is the magnitude of vector xx.
%
% CALLED FUNCTIONS:
%   none
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
%   01/01/2017, GIANMARIO MERISIO, Added case i=0, i=pi, e=0, e=1;
%
% -------------------------------------------------------------------------

kep=zeros(1,6);
rr=state(1:3);
vv=state(4:6);
r = norm(rr,2);
v = norm(vv,2);

TOL=1e-9; %value of tolerance due to ode options

% Compute h
hh = cross(rr,vv);
h = norm(hh,2);

% Compute e
ee = 1/mu * cross(vv,hh) - rr./r;
e = norm(ee,2);
ee_v = ee./e;       %versor

if abs(1-e)<TOL           % Parabolic orbit
    kep(1)=h;
else
    % Compute a
    a = (2./r - v.*v/mu).^(-1);
    kep(1)=a;
end

% Compute i
i = acos(hh(3)./h);
kep(3)=i;

% Compute Om 
if i==0 || i==pi   % Equatorial orbit, Pro/Retrograde
    Om=0;
    NN=[1 0 0];
else               % Inclined orbit
    NN = cross([0 0 1],hh);
    N= norm(NN,2);
    NN= NN./N;
    Om = acos(NN(1));
    if(NN(2)<0)
        Om = 2*pi - Om;
    end
end
kep(4)=Om;


% Compute w
if i==0 || i==pi        % Equatorial orbit, Pro/Retrograde
    if e < TOL          % Circular Orbit
        w=0;
    else                % NON-Circular Orbit            
        w = acos(ee_v(1));
        if i==0 && ee_v(2)<0                
            w = 2*pi - w;
        elseif i==pi && ee_v(2)>0               
            w = 2*pi - w;
        end
    end
else                    % Inclined orbit
    if e < TOL          % Circular Orbit
        w=0;
    else                % NON-Circular Orbit            
        w = acos(NN*ee_v');
        if ee_v(3)<0
            w = 2*pi - w;
        end
    end
end
kep(5)=w;


% Compute theta
if i==0 || i==pi        % Equatorial orbit, Pro/Retrograde
    if e < TOL          % Circular Orbit
        theta = acos(rr(1)/r);          
        if i==0 && rr(2) < 0                    
            theta = 2*pi - theta;
        elseif i==pi && rr(2) > 0                    
            theta = 2*pi - theta;
        end   
    else                % NON-Circular Orbit            
        re = rr*ee_v';
        frac=re/r;
        if frac<-1      % To avoid complex value,sometimes they appear
            frac=-1;
        elseif frac>1
            frac=1;
        end
        theta = acos(frac);
        if vv*rr'<0
            theta = 2*pi - theta;
        end
    end
else                    % Inclined orbit
    if e < TOL          % Circular Orbit
        rn = rr*NN';
        frac=rn/r;
        if frac<-1      % To avoid complex value,sometimes they appear
            frac=-1;
        elseif frac>1
            frac=1;
        end
        theta = acos(frac);
        if rr(3)<0  
            theta = 2*pi - theta;
        end   
    else                % NON-Circular Orbit            
        re = rr*ee_v';
        frac=re/r;
        if frac<-1      % To avoid complex value,sometimes they appear
            frac=-1;
        elseif frac>1
            frac=1;
        end
        theta = acos(frac);
        if vv*rr'<0
            theta = 2*pi - theta; 
        end
    end
end
kep(6)=theta;

if e<TOL
    e=0;                % Circular orbit
elseif abs(1-e)<TOL
    e=1;                % Parabolic orbit
end
kep(2)=e; % This is the last thing because I adjust eccentricity due to different cases.

return