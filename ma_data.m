planet_1 = 3;                           % Earth number according to uplanet function
planet_2 = 4;                           % Mars number according to uplanet function
mu_earth = astroConstants(13);          % Earth's planetary constant 
mu_mars = astroConstants(14);           % Mars' planetary constant
mu_sun = astroConstants(4);             % Sun planetary constant
mu_phobos = 0.0007112;                  % Phobos gravitational constant
G = astroConstants(1);                  % Universal gravitational constant
AU = astroConstants(2);                 % Astronomic Unit
m_mars = mu_mars/G;                     % Mass of Mars
m_sun = mu_sun/G;                       % Mass of the Sun
m_earth = mu_earth/G;                   % Mass of the Earth
m_phobos = mu_phobos/G;                 % Mass of Phobos
mr_earth = astroConstants(23);          % Earth's mean radius
mr_mars = astroConstants(24);           % Mars' mean radius
options = odeset('Reltol',1e-13,'AbsTol',1e-14);