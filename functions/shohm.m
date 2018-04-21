function    [DV1, DV2, DV_TOT] = shohm(kep_1,kep_2,mu)

% shohm.m -    Returns departure, arrival and total DVs in case of a
%              simil-Hohmann approximantion of the transfer arc.
%
% PROTOTYPE:
%   [DV1, DV2, DV_TOT] = shohm(kep_1,kep_2,mu)
%
% DESCRIPTION:
%              This function approximate the orbits as circular ones.
%              Radius of each orbit is choosen as the modulus of the 
%              radius vector in the heliocentric inertial frame.
%              Orbits are supposed to be coplanar.
%
%   List of identifiers:
%       
%
%   Notes for upgrading this function:
%       
%
% INPUT:    
%             -kep_1 -> keplerian parameters of first celestial body 
%             -kep_2 -> keplerian parameters of second celestial body
%             -mu    -> gravitational constant of the main attractor
%
% OUTPUT:   
%             -DV1    -> departure DV
%             -DV2    -> arrival DV
%             -DV_TOT -> total DV 
%   
% EXAMPLE:
%
% REFERENCES:
%
% CALLED FUNCTIONS:
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
% PREVIOUS VERSION:
%   29/12/2016, ALESSANDRO MARIA MASSERINI
%
% CHANGELOG:
%
% -------------------------------------------------------------------------
state_1=kep2car(kep_1,mu);
state_2=kep2car(kep_2,mu);

rr_1=kep_1(1:3);
rr_2=kep_2(1:3);
r_1=norm(rr_1);
r_2=norm(rr_2);

rp=r_1;
ra=r_1+r_2;
a=(rp+ra)/2;
e=(r_2-r_1)/(r_2+r_1);
P=a*(1-e^2);

vp=sqrt(mu/P)*(1+e);
va=sqrt(mu/P)*(1-e);
v_1=sqrt(mu/r_1);
v_2=sqrt(mu/r_2);

DV1=abs(vp-v_1);
DV2=abs(v_2-va);
DV_TOT=DV1+DV2;
