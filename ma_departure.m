%% Definition of the launch and arrival window----------------------------- 
date_depi = [2023 07 01 12 00 00];      % Initial departure date
date_depf = [2025 09 30 12 00 00];      % Final departure date
date_arri = [2024 05 01 12 00 00];      % Initial arrival date
date_arrf = [2026 08 31 12 00 00];      % Final arrival date

% Conversion to mjd2000----------------------------------------------------
date_depi_mjd2000 = date2mjd2000(date_depi);    % Initial departure date [mjd2000]
date_depf_mjd2000 = date2mjd2000(date_depf);    % Final departure date   [mjd2000]
date_arri_mjd2000 = date2mjd2000(date_arri);    % Initial arrival date   [mjd2000]
date_arrf_mjd2000 = date2mjd2000(date_arrf);    % Final arrival date     [mjd2000]

% Definition of launch and arrival window arrays---------------------------
N = 5;                                              % Days per step 
dep_dates = date_depi_mjd2000:N:date_depf_mjd2000;  % Departure dates vector
arr_dates = date_arri_mjd2000:N:date_arrf_mjd2000;  % Arrival dates vector

%% Pork chop plot of the Lambert's problem----------------------------------
[PCdata] = porkChopDatav3(planet_1,planet_2,dep_dates,arr_dates);
save PCdata.mat
dep=PCdata.dep;
arr=PCdata.arr;
TOF=PCdata.TOF;
dv1=PCdata.dv1;
dv2=PCdata.dv2;
DV=PCdata.DV;
dv1_dMAX=3.5;
porkChopPlot(PCdata,1,0,1,dv1_dMAX);
