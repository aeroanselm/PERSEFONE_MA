function [hfig] = drawPlanet( namePlanet, position, handle, scale)
% drawPlanet.m - draws a planet surface (approx) in a given position.
%
% PROTOTYPE: 
%           [HFIG] = drawPlanet( namePlanet, position, varargin)
%
% DESCRIPTION: 
%  The function draws a planet in a given position and also a
%  given axis.
%
% INPUT:
%  namePlanet     -    string with the name of the planet (not case
%                           sensitive)
%  position       -    centre of the planet to draw
%  handle         -    Figurr or axis handle
%
% OUTPUT:
%  HFIG           -    Figure object of the planet draw.
%
% EXAMPLE: hfig = drawPlanet('Sun',[0 0 0]);
%          hfig = drawPlanet('Earth',[0,0,0])
%          Note: you need to add the folder textures to the Matlab path!
%
% CALLED FUNCTIONS: getAstroConstants.m, folder textures
%
% AUTHOR:
%   Joan Pau Sanchez, 14/10/2009, MATLAB, drawPlanet.m
%   
% CHANGELOG:
%   Camilla Colombo, 31/10/2016: example added, tidied up.
%   REVISION: 31/10/2016, Camilla Colombo
%   Gianmario Merisio, 03/11/2016: now it works for current folder, change
%                                  getAstroConstants with astroConstants,
%                                  it does not work with 'Asteroid'
%
%--------------------------------------------------------------------------

% 1. Searching Texture Directory
persistent pathTextures
sep = filesep;
if isempty(pathTextures)
    %p = path;
    p=[pwd sep 'functions' sep 'texture'];
    remain=p;
    token='a';
    i=1;
    pathTextures=[];
    while isempty(pathTextures) && strcmp(token,'')~=1
        [token, remain] = strtok(remain, ';');
        %x = strfind(token, 'textures');
        x = strfind(token, 'texture');
        if any(x)
            pathTextures=token;
        end
    end
    if isempty(pathTextures)
        msg = ['The texture files were not found in the path'];
        eid = sprintf('TOOLBOX:%s:propertyError', mfilename);
        error(eid,'%s',msg)
    end
end
%--------------------------------------------------------------------------
% 2. Loading Texture Map 
     files = dir(pathTextures);
     index = strfind({files(:).name}, namePlanet);
    texturemap = [];
    i=1;
    while isempty(texturemap)
        if isempty(index{i})
            i=i+1;  
        else
            [texturemap,map] = imread([pathTextures sep files(i).name]);
        end
    end
%--------------------------------------------------------------------------
% 3. Preraring Figure Object
if nargin<3
    HAXIS = gca;
elseif ishandle(handle)==0
        msg = ['The figure handle is not valid'];
        eid = sprintf('TOOLBOX:%s:propertyError', mfilename);
        error(eid,'%s',msg)
else
    try
        HAXIS=gca(handle);
    catch
        HAXIS=handle;  
    end
    hold on
end
%--------------------------------------------------------------------------
if nargin<4
    scale=1;   
end
% 4. Plotting Planet
switch namePlanet
    case 'Mercury'
        numP=21;
    case 'Venus'
        numP=22;
    case 'Earth'
        numP=23;
    case 'Mars'
        numP=24;
    case 'Jupiter'
        numP=25;
    case 'Saturn'
        numP=26;
    case 'Uranus'
        numP=27;
    case 'Neptune'
        numP=28;
    case 'Moon'
        numP=30;
    case 'Sun'
        numP=3;
    case 'Asteroid'
        numP=31;
    case 'Phobos'
        numP=32;
    otherwise
        error('This planet is missing.');
end
    
%[radPlanet] = getAstroConstants(namePlanet, 'Radius');
if numP~=31 
    [radPlanet] = astroConstants(numP);
else
    radPlanet=0.6;      % Radius of Hermes, note that is a binary asteroid
end
    nfaces = 50;
    [X,Y,Z] = sphere(nfaces);
% 4.1. scaling the sphere and locating the planet
X = -radPlanet*X*scale + position(1);
Y = -radPlanet*Y*scale + position(2);
Z = -radPlanet*Z*scale + position(3);
% 4.2. plotting
surf(HAXIS, X, Y, Z,texturemap, 'LineStyle', 'none','FaceColor', 'texturemap');
hfig = gcf; 
return                     