function x_d = dyn_2BP (t, x, mu)

% dyn_2BP.m - Dynamics of 2 Body Problem.
% 
% TYPE: 
%   Function.
%
% PROTOTYPE:
%   x_d = dyn_2BP (t, x, mu)
%
% DESCRIPTION:
%   Function that contains the dynamics of the 2 Body Problem. It is used
%   when you call ode45 or ode113. It is preferable to use it with ode113
%   when you are working with astrodynamics.
%
%   Note:   Use it with proper options. e.g.:
%           options=odeset('RelTol',1e-13,'AbsTol',1e-14);
%           [t,state]=ode113(@(t,y)dyn_2BP(t,y,mu), time, state0, options);
%   
%   See also:   plotOrbit.m
%               ode45.m
%               ode113.m
%
% INPUT:
%	t[1]    Time instant.    
%   x[6]    State variables of the system.
%   mu[1]   Planetary constant of the planet (mu = mass * G) [L^3/T^2]
%   
% OUTPUT:
%   x_d[6]    Derivatives of the state.
%
% NOTATION:
%   kep[6] is a vector made as kep=[a,e,i,Om,w,theta].
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
%   dd/mm/yyyy, NAME SURNAME, informations
%
% -------------------------------------------------------------------------

[n,m]=size(x);
x_d=zeros(n,m);
r=sqrt(x(1)^2+x(2)^2+x(3)^2);
a=-mu/r^3;

x_d(1)=x(4);
x_d(2)=x(5);
x_d(3)=x(6);
x_d(4)=a*x(1);
x_d(5)=a*x(2);
x_d(6)=a*x(3);

end


