function [dv] = ffdv (X,flag,num1,num2)

% ffdv.m - Fitness Function of dv at Departure, Arrival or total.
% 
% TYPE: 
%   Function.
%
% PROTOTYPE:
%   [dv] = ffdv (X,flag,num1,num2)
%
% DESCRIPTION:
%   This is the fitness function for the computation of dv at departure,
%   arrival or total, it depends on the value of flag. This function has been
%   written with the purpose to be used wiht the function 'ga'. 'ga' is a function
%   that applies genetic algorithm to find a minimum. This function is
%   specific for interplanetary mission in the Solar System, so mu is
%   assumed to be the mu of the Sun!
%
%   Note: 
%   
%   See also:   ga.m
%               lambertMR.m
%
% INPUT:
%	X[2]    Variables of the problem, the are 2 dates in the format of mjd2000.
%           The first one is the departure date, the last one is the
%           arrival date.
%   flag[1] This flag is used to choose which dv value has to be returned:
%               0 -> dv at departure
%               1 -> dv at arrival
%               2 -> dv total
%   num1[1] Number of the celestial body 1.
%   num2[1] Number of the celestial body 2.
%
% OUTPUT:
%   dv[1]    Scalar return of the fitness function. It is dv at arrival or
%            departure depending on the flag value [Km/s].
%
% NOTATION:
%   kep[6] is a vector made as kep=[a,e,i,Om,w,theta].
%   state[6] is a vector made as state=[x,y,z,vx,vy,vz].
%   xx is a vector.
%   x is the magnitude of vector xx.
%
% CALLED FUNCTIONS:
%   lambertMR.m
%   ephNEO.m
%   uplanet.m
%   astroConstants.m
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

day1=X(1);
day2=X(2);
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
%dv_tot=norm(state1(4:6)-state1_tran(4:6))+norm(state2(4:6)-state2_tran(4:6));

if flag==0
    dv=(norm(state1_tran(4:6)-state1(4:6)));
elseif flag==1
    dv=(norm(state2_tran(4:6)-state2(4:6)));
else
    dv=(norm(state1_tran(4:6)-state1(4:6)))+(norm(state2_tran(4:6)-state2(4:6)));
end

    
return


