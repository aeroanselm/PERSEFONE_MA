function [jd] = mjd20002jd (mjd2000)

% mjd20002jd.m - Julian day number from Modified Julian day 2000 number to.
% 
% TYPE: 
%   Function.
%
% PROTOTYPE:
%   [jd] = mjd20002jd (mjd2000)
%
% DESCRIPTION:
%   It returns the Julian day number from the Modified Julian day 2000
%   number.
%
%   Note: 
%   
%   See also: 
%
% INPUT:
%   mjd2000[1]  Date in MJD 2000. MJD2000 is defined as the number of days
%               since 01-01-2000, 12:00 noon. It must be a real greater or
%               equal than -2451545.     
%
% OUTPUT:
%   jd[1]       Date in Julian Day. The JD (Julian day) count is from 0 at
%               12:00 noon, 1 January -4712 (4713 BC), Julian proleptic
%               calendar. The corresponding date in Gregorian calendar is
%               12:00 noon, 24 November -4713.
%
% NOTATION:
%
% CALLED FUNCTIONS:
%   none
%
%
% ORIGINAL VERSION:
%   16/02/2008, Nicola Croisard, this function was missing from the toolbox,
%   so I rewrite it.
%
% AUTHOR
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

jd = mjd2000 + 2400000.5 + 51544.5;

return


