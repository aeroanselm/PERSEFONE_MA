function [dv_d,dv_a] = evaldv (DAY,num1,num2)

% evaldv.m - Evaluation of dv at departure and arrival.
% 
% TYPE: 
%   Function.
%
% PROTOTYPE:
%   [dv_d,dv_a] = evaldv (day1,day2,num1,num2)
%
% DESCRIPTION:
%   This function evaluates dv_d and dv_a for a given interplanetary mission.
%   It is supposed that the interplanetary mission is inside the Solar System.
%
%   Note:   it is very usefull call this function after you used algorithms
%           that give you a minimum for the lambert problem.
%   
%   See also: lambertMR.m
%
% INPUT:
%	DAY[2]  2 dates in the format of mjd2000. The first one is the departure 
%           date, the last one is the arrival date.
%   num1[1] Number of the celestial body 1.
%   num2[1] Number of the celestial body 2.     
%
% OUTPUT:
%   dv_d[1]  Value of dv at the departure [Km/s].
%   dv_a[1]  Value of dv at the arrival [Km/s].
%
% NOTATION:
%   kep[6] is a vector made as kep=[a,e,i,Om,w,theta].
%   state[6] is a vector made as state=[x,y,z,vx,vy,vz].
%   xx is a vector.
%   x is the magnitude of vector xx.
%
% CALLED FUNCTIONS:
%   stroConstants.m
%   ephNEO.m
%   lambertMR.m
%   kep2car.m
%
% REFERENCES:
%   - Author,"Book's Name or Article", dd/mm/yyyy.
%
% FUTURE DEVELOPMENT:
%   Work in progress.
%
% ORIGINAL VERSION:
%   dd/mm/yyyy, NAME SURNAME, informations.
%
%--------------------------------------------------------------------------
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

% mu Sun
mu=astroConstants(4);

day1=DAY(1);
day2=DAY(2);
dof=day2-day1;
tof=dof*24*3600;
if num1>=12
    kep1=ephNEO(day1,num1);
else
    kep1=uplanet(day1,num1);
end
if num2>=12
    kep2=ephNEO(day2,num2);
else
    kep2=uplanet(day2,num2);
end
state1=kep2car(kep1,mu);
state2=kep2car(kep2,mu);
[~,~,~,~,VI,VF,~,~] = lambertMR(state1(1:3),state2(1:3),tof,mu,0,0,0);
state1_tran=[state1(1:3) VI];
state2_tran=[state2(1:3) VF];

dv_d=(norm(state1_tran(4:6)-state1(4:6)));
dv_a=(norm(state2_tran(4:6)-state2(4:6)));


return


