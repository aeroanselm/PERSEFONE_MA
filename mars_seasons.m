clear 
close all
clc

addpath(genpath('./functions'))         % Adding the functions folder to path [linux version, works even on Windows]
date_i = [2026 01 01 12 00 00];
date_i = date2mjd2000(date_i);
year = 365;
days = date_i:(date_i + 3*year);

figure ()
hold on
grid on

for i = days
    kepseason = uplanet(i,4);
    [season,Ls] = mseasons(kepseason(6)*180/pi);
    switch season
        case 'Spring'
            plot(i,1,'oy');
        case 'Summer'
            plot(i,1,'or');
        case 'Fall'
            plot(i,1,'og');
        case 'Winter'
            plot(i,1,'ob');
    end
end