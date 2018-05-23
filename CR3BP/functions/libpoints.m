% apo_peri.m - Returns departure, arrival and total DVs in case of a
%              bitangent approximantion of the transfer arc.
%
% PROTOTYPE:
%   [ap_DV1, ap_DV2, ap_DV_TOT]=apo_peri(kep_1,kep_2,mu)
%
% DESCRIPTION:
%              In this case the departure point is set on the apoapsis of the 
%              first orbit and the arrival point on the periapsis of the
%              second  one.
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
%             -ap_DV1    -> departure DV
%             -ap_DV2    -> arrival DV
%             -ap_DV_TOT -> total DV 
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

function [L,d] = libpoints(m1,m2,r12)

M = m1+m2;
G = astroConstants(1);
mu = G*M;
p1 = m1/(m1+m2);
p2 = m2/(m1+m2);

% Synodic reference
d1 = -p2*r12;
d2 = p1*r12;
d = [d1; d2];

fzeta = @(zeta) ((1-p2)./abs(zeta+p2).^3) .* (zeta+p2) + p2./abs(zeta+p2-1).^3 .* (zeta+p2-1)-zeta;

x1 = fzero(fzeta,0.8)*r12;  y1 = 0;
x2 = fzero(fzeta,1.1)*r12;  y2 = 0;
x3 = fzero(fzeta,-1.1)*r12; y3 = 0;
x4 = r12/2-p2*r12;          y4 = +sqrt(3)/2*r12;
x5 = x4;                    y5 = -sqrt(3)/2*r12;

L = [x1 y1; x2 y2; x3 y3; x4 y4; x5 y5];
end

