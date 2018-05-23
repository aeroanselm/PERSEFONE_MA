function [] = PCplot (dep,arr,tof,A,flag,name,uom,ulim)

% PCplot.m - Pork Chop Plot from a set of data.
% 
% TYPE: 
%   Function.
%
% PROTOTYPE:
%   [] = PCplot (dep,arr,tof,A,flag,name,uom,ulim)
%
% DESCRIPTION:
%   This Routine plots the pork chop plot from a set of data. The plot can
%   be use to analyse and to search suitable solutions for interplanetary
%   missions.
%
%   Note: 
%   
%   See also:   porkChopData.m
%               PCplot.m
%
% INPUT:
%   dep[y]          Starting days for the transfer [MJD2000].
%   arr[x]          Arrival days for the transfer [MJD2000].
%   tof[x,y]        Times of flight [days].
%   A[x,y]          Generic quantity [-].
%   flag[1]         It determines different type of axes:
%                       0 ->    MJD2000
%                       1 ->    Classical Date
%   name[string]    Name of A.
%   uom[string]     Unit of measure of A.
%   ulim[1]         Upper limit for level lines (OPTIONAL).
%
% OUTPUT:
%   []
%
% NOTATION:
%   kep[6] is a vector made as kep=[a,e,i,Om,w,theta].
%   state[6] is a vector made as state=[x,y,z,vx,vy,vz].
%   xx is a vector.
%   x is the magnitude of vector xx.
%
% CALLED FUNCTIONS:
%   mjd20002date.m
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

if nargin==7
    ulim=60;                                    %SET
end
pitch=1;                                   %SET
Ls=1000;                                    %SET

switch flag
    case 0
    case 1
        DEP=zeros(length(dep),6);
        ARR=zeros(length(arr),6);
        for it=1:length(dep)
            DEP(it,:)=mjd20002date(dep(it));
        end
        for it=1:length(arr)
            ARR(it,:)=mjd20002date(arr(it));
        end
        clear dep;
        clear arr;
        dep=datenum(DEP);
        arr=datenum(ARR);
    otherwise
        error ('WRONG flag.');
end

figure();
hold on;
% A curves
Alevs=0:pitch:ulim;                                     %SET
Alevs(1)=1;
contour(dep,arr,A,Alevs);
Aindex=0:10:ulim;
Aindex(1)=1;
[C,h]=contour(dep,arr,A,Aindex,'lineWidth',2);
clabel(C,h,Aindex(1:4),'FontSize',15);
colorbar;
caxis([0,ulim]);
        
% Time of Flight lines
toflevs=0:250:max(max(tof));                           %SET
toflevs(1)=10;
[C,h]=contour(dep,arr,tof,toflevs,'r','LineWidth',1.5);
clabel(C,h,'FontSize',15,'Color','r','LabelSpacing',Ls);

legend([name ' ' uom],[name ' ' uom],'Time of Flight [days]','Location','southeast');


switch flag
    case 0
        xlabel(['DEPARTURE [MJD2000]']);
        xtickangle(45);
        xticks(dep(1):30:dep(end));
        ylabel(['ARRIVAL [MJD2000]']);
        yticks(arr(1):50:arr(end)); 
        title(['Pork Chop Plot, ' name ' ' uom]);
    case 1
%         datetick('y',29,'keeplimits','keepticks');
%         datetick('x',29,'keeplimits','keepticks');
        xticks(dep(1):30:dep(end));                               %SET
        xtickangle(45);
        yticks(arr(1):50:arr(end));                               %SET
        datetick('y',1,'keeplimits','keepticks');
        datetick('x',1,'keeplimits','keepticks');
        xlabel(['DEPARTURE']);
        ylabel(['ARRIVAL']);
        title(['Pork Chop Plot, ' name ' ' uom]);
    otherwise
        error ('WRONG flag.');
end
hold off;

return


