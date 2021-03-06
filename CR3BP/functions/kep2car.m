function [state]= kep2car (kep,mu)

% kep2car.m - From keplerian to cartesian coordinates.
% 
% TYPE: 
%   Function.
%
% PROTOTYPE:
%   [state]= kep2car (kep,mu)
%
% DESCRIPTION:
%   The function uses rotational matrices to transform the keplerian values
%   of the orbit into the cartesian coordinates in the inertial frame.
%   
%   See also:   car2kep.m
%
% INPUT:
%   kep[6]      Keplerian values of the orbit.
%   mu[1]       Planetary constant of the planet (mu = mass * G) [L^3/T^2]          
%
% OUTPUT:
%   state[6]    It is the state of the spacecraft in the cartesian
%               reference frame
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
%   01/01/2017, GIANMARIO MERISIO, Added case e=1.
%
% -------------------------------------------------------------------------


e=kep(2);
if e==1
    h=kep(1);
else
    a=kep(1);
end
i=kep(3);
Om=kep(4);
w=kep(5);
theta=kep(6);

% ROTATIONAL MATRIX

%[Om]_k
R3 = [cos(Om), sin(Om), 0;
     -sin(Om), cos(Om), 0;
      0,      0,        1];
%[i]_i
R2 = [1     ,0        ,0;
      0     ,cos(i)   ,sin(i);
      0     ,-sin(i)  ,cos(i)];
%[w]_k
R1 = [cos(w)  ,sin(w) ,0;
      -sin(w) ,cos(w) ,0;
       0      ,0      ,1];

R = [R1*R2*R3]';

% compute distances rho = rho(theta) fixed a,e or h,e (parabola)
if e==1
    p=h^2/mu;
else
    p = a*(1-e*e);
end
rho = p./(1+e*cos(theta));

%collect r_PE,v_PE 
r_PE = [rho.*cos(theta)  rho.*sin(theta)  0];
v_PE = (mu/p)^(0.5).*[-sin(theta)  e+cos(theta) 0];

state = [(R*r_PE')' (R*v_PE')'];

return
