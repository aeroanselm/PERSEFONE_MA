function [] = porkChopPlot(porkChop_data,option,flag,flagcon,constraint,interest,ulim)	

% porkChopPlot.m - Pork Chop Plot from a struct of data.
%  
% TYPE: 
%   Function.
%
% PROTOTYPE:
%   [] = porkChopPlot(porkChop_data,option,flag,flagcon,constrain,ulim)
%
% DESCRIPTION:
%   This Routine plots the pork chop plot from a set of data. The plot can
%   be use to analyse and to search suitable solutions for interplanetary
%   missions. This function was written with the idea of a constraint on
%   the C3d. Of course you can eventually change the function so that it
%   will work with a different type of constraint. In the places where
%   there is the world 'SET' you can choose your own values. Specific
%   routine for the first assignment.
%
%   Note: 
%   
%   See also:   porkChopData.m
%
% INPUT:
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
%   option[1]               It determines different types of the plot
%                           content:
%                               0 ->    C3 Pork Chop Plot
%                               1 ->    dv Pork Chop Plot 
%   flag[1]                 It determines different types of axes:
%                               0 ->    MJD2000
%                               1 ->    Classical Date
%   flagcon[1]              It determines the type of the constrain:
%                               1 ->    on departure value (C31/dv1)
%                               2 ->    on arrival value (C32/dv2)
%                               3 ->    on total value (C3/DV)
%   constraint[1]           Define a constraint of the problem.
%   interest[1]             Activate specific plot for Ass.1 [optional] :
%                               0 ->    off
%                               1 ->    on 
%   ulim[1]                 Upper limit for level lines (OPTIONAL).
%
% OUTPUT:
%   []
%   The output is the Pork Chop plot.
%
% NOTATION:
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
%   ALESSANDRO MARIA MASSERINI, creation of this function - 17/12/2016
%   GIANMARIO MERISIO, change of plot style - 19/12/2016
% --------------------------------------------------------------------------

if nargin <=5
    interest=0;
end

dep=porkChop_data.dep;
arr=porkChop_data.arr;
if option==0
    Qtot=porkChop_data.C3;
    Q1=porkChop_data.C31;
    Q2=porkChop_data.C32;
else 
    Qtot=porkChop_data.DV;
    Q1=porkChop_data.dv1;
    Q2=porkChop_data.dv2;
end
TOF=porkChop_data.TOF;

Ls=1000;                                %SET

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
        dep=datenum(DEP);
        arr=datenum(ARR);
    otherwise
        error ('WRONG flag.');
end

figure();
hold on;
if interest==0
    % Q1 curves
    Q1levs=0:2:30;                                     %SET
    Q1levs(1)=1;
    contour(dep,arr,Q1,Q1levs,'g');

    % Q2 curves
    Q2levs=0:5:30;                                     %SET
    Q2levs(1)=1;
    contour(dep,arr,Q2,Q2levs,'y');

    % C3 curves
    Qtotlevs=0:2:30;                                     %SET
    Qtotlevs(1)=1;
    contour(dep,arr,Qtot,Qtotlevs,'k');
else
    if nargin == 6
        ulim=60;                                                    %SET
    end
    pitch=1;                                                    %SET
    % Qtot curves
    Qtotlevs=0:pitch:ulim;                                     %SET
    Qtotlevs(1)=1;
    contour(dep,arr,Qtot,Qtotlevs);
    Qtotindex=0:10:ulim;
    Qtotindex(1)=1;
    [C,h]=contour(dep,arr,Qtot,Qtotindex,'lineWidth',2);
    clabel(C,h,Qtotindex(1:4),'FontSize',15);
    colorbar;
    caxis([0,ulim]);
end
% Constraint
switch flagcon
    case 1
        [C,h]=contour(dep,arr,Q1,constraint*[1 1],'r','LineWidth',2);
        clabel(C,h,'FontSize',12,'LabelSpacing',Ls);
    case 2
        [C,h]=contour(dep,arr,Q2,constraint*[1 1],'r','LineWidth',2);
        clabel(C,h,'FontSize',12,'LabelSpacing',Ls);
    case 3
        [C,h]=contour(dep,arr,Qtot,constraint*[1 1],'r','LineWidth',2);
        clabel(C,h,'FontSize',12,'LabelSpacing',Ls);
    otherwise
        error('flagcon value is wrong!');
end

% Time of Flight lines
toflevs=0:250:max(max(TOF));                    %SET
toflevs(1)=10;
[C,h]=contour(dep,arr,TOF,toflevs,'b','LineWidth',2);
clabel(C,h,'FontSize',12,'Color','b','LabelSpacing',Ls);

if option==0
    if interest==0
        legend('C3 departure [Km^2/s^2]', 'C3 arrival [Km^2/s^2]', 'C3_{tot} [Km^2/s^2]', 'Constraint [Km^2/s^2]', ...
            'Time of Flight [days]','Location','southeast');
    else
        legend('C3_{tot} [Km^2/s^2]','C3_{tot} [Km^2/s^2]','Constraint [Km^2/s^2]','Time of Flight [days]','Location','southeast');
    end
else
    if interest==0
        legend('\Deltav departure [Km/s]', '\Deltav arrival [Km/s]', '\DeltaV_{tot} [Km/s]', 'Constraint [Km/s]', ...
            'Time of Flight [days]','Location','southeast');
    else
        legend('\DeltaV_{tot} [Km/s]','\DeltaV_{tot} [Km/s]','Constraint [Km/s]','Time of Flight [days]','Location','southeast');
    end
end
switch flag
    case 0
        xlabel('DEPARTURE [MJD2000]','Interpreter','latex');
        xtickangle(45);
        xticks(dep(1):50:dep(end));
        ylabel('ARRIVAL [MJD2000]','Interpreter','latex');
        yticks(arr(1):50:arr(end)); 
        title('\textbf{Pork Chop Plot, destination: Earth}','Interpreter','latex');
    case 1
%         datetick('y',29,'keeplimits','keepticks');
%         datetick('x',29,'keeplimits','keepticks');
        xticks(dep(1):30:dep(end));                               %SET
        xtickangle(45);
        yticks(arr(1):50:arr(end));                               %SET
        datetick('y',1,'keeplimits','keepticks');
        datetick('x',1,'keeplimits','keepticks');
        xlabel('DEPARTURE','Interpreter','latex');
        ylabel('ARRIVAL','Interpreter','latex');
        title('\textbf{Pork Chop Plot, destination: Earth}','Interpreter','latex');
    otherwise
        error ('WRONG flag.');
end
hold off;

return





