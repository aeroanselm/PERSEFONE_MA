function [porkChop_data] = porkChopData(p1,p2,date_li,date_lf,date_ai,date_af,N)

% porkChopData.m - Data needed for Pork Chop plot.
% 
% TYPE: 
%   Function.
%
% PROTOTYPE:
%   [porkChop_data] = porkChopData(p1,p2,date_li,date_lf,date_ai,date_af,N)
%
% DESCRIPTION:
%   The function computes the value of several quantities for each combination 
%   of the interval of date and time of flight. This function  uses the 
%   ephemerides. It is a specific function for interplanetary missions.
%   Some values are specific of the solar system and the Sun.
%   
%   See also:   lambertMR.m
%               porkChopPlot.m
%
% INPUT:
%   p1[1]       Number of the celestial body 1.
%   p2[1]       Number of the celestial body 2.
%   date_li[6]  INITIAL date, it is referring to the first point.
%               Date in the Gregorian calendar, as a 6-element vector
%               [year, month, day, hour, minute, second]. For dates before
%               1582, the resulting date components are valid only in the
%               Gregorian proleptic calendar. This is based on the
%               Gregorian calendar but extended to cover dates before its
%               introduction. date must be after 12:00 noon, 24 November
%               -4713.
%   date_lf[6]  FINAL date, it is referring to the first point defined as date_li.
%   date_ai[6]  INITIAL date, it is referring to the second point defined as date_li.
%   date_af[6]  FINAL date, it is referring to the second point defined as date_li.
%   N[1]        Pitch for the matrices of departure window and arrival window.
%
% OUTPUT:
%   porkChop_data [struct]  It the result of the calculation of the Lambert's 
%                           problem for the launch window. It contains:
%                               dep[y]      Starting days for the transfer [MJD2000].
%                               arr[x]      Arrival days for the transfer [MJD2000].
%                               TOF[x,y]    Times of flight [days].
%                               dv1[x,y]    Cost of the manoeuvre at departure [Km^2/s^2].
%                               dv2[x,y]    Cost of the manoeuvre at arrival [Km^2/s^2].
%                               DV[x,y]     Total Cost of the manoeuvres [Km^2/s^2].
%                               C3[x,y]     Characteristic energy, C3d+C3a [Km^2/s^2].
%                               C31[x,y]    Characteristic energy for interplanetary 
%                                           missions at departure [Km^2/s^2].
%                               C32[x,y]    Characteristic energy for interplanetary 
%                                           missions at arrival [Km^2/s^2].
%
% NOTATION:
%   kep[6] is a vector made as kep=[a,e,i,Om,w,theta].
%   state[6] is a vector made as state=[x,y,z,vx,vy,vz].
%   xx is a vector.
%   x is the magnitude of vector xx.
%
% CALLED FUNCTIONS:
%   date2mjd2000.m
%   kep2car.m
%   lambertMR.m
%   uplanet.m
%   ephNEO.m
%   astroConstants.m
%
% FUTURE DEVELOPMENT:
%   Work in progress.
%
% REFERENCES:
%   - Author,"Book's Name or Article", dd/mm/yyyy.
%
% FUTURE DEVELOPMENT:
%   Work in progress.
%
% ORIGINAL VERSION:
%   17/12/2016, ALESSANDRO MARIA MASSERINI
%				GIANMARIO MERISIO
%				GIOVANNI ZANOTTI
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
%   ALESSANDRO MARIA MASSERINI creation of this function - 17/12/2016
%
% --------------------------------------------------------------------------

mu_s=astroConstants(4);				% Sun gravitational constant
day=24*3600;						% Conversion from days to seconds

% Conversion from Gregorian dates to mjd2000
mjd2000_li=date2mjd2000(date_li);
mjd2000_lf=date2mjd2000(date_lf);
mjd2000_ai=date2mjd2000(date_ai);
mjd2000_af=date2mjd2000(date_af);

dep=[mjd2000_li:N:mjd2000_lf];		% Launch date vector
arr=[mjd2000_ai:N:mjd2000_af];		% Arrival date vector

DV=zeros(length(dep),length(arr));	% DV matrix initialisation
dv1=DV; dv2=DV; TOF=DV;				% Other output matrices initialisation	

for i=1:length(dep)
	for j=1:length(arr)
		dep_date=dep(i);       		% Departure date
        	arr_date=arr(j);        % Arrival date
		
		if dep_date<arr_date
            
            if p1<11
                kep_p1=uplanet(dep_date,p1);   	% Ephemerides call for p1
            else
                kep_p1=ephNEO(dep_date,p1);
            end
            
			if p2<11
				kep_p2=uplanet(arr_date,p2);	% Ephemerides call for p2
			else
				kep_p2=ephNEO(arr_date,p2);		% Ephemerides call for p2
			end
			
			car_p1=kep2car(kep_p1,mu_s);        % State vector for p1
    		car_p2=kep2car(kep_p2,mu_s);		% State vector for p2
        		
        	R1=car_p1(1:3);						% Position vector for p1
        	R2=car_p2(1:3);						% Position vector for p2
        		
        	tof=(arr_date-dep_date);			% Time of flight expressed in days
        	tof_sec=(arr_date-dep_date)*day;    % Time of flight expressed in seconds
        		
        	[~,~,~,~,V1,V2]=lambertMR(R1,R2,tof_sec,mu_s,0,0,0,0);	% Solving Lambert's problem
        		
        	dv1(i,j)=norm(V1-car_p1(4:6));
        	dv2(i,j)=norm(car_p2(4:6)-V2);
			TOF(i,j)=tof;
		else
        	dv1(i,j)=NaN;
        	dv2(i,j)=NaN;
			TOF(i,j)=NaN;	
		end
	end	
end

dv1=dv1';
dv2=dv2';
TOF=TOF';
DV=dv1+dv2;
C31=dv1.^2;
C32=dv2.^2;
C3=C31+C32;

porkChop_data=struct('dep',dep,'arr',arr,'TOF',TOF,'dv1',dv1,'dv2',dv2,...
                    'DV',DV,'C31',C31,'C32',C32,'C3',C3);

return
		

