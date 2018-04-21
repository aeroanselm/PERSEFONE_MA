EC= 1.425189346826480E-01;
OM= 2.356326194919086E+02*pi/180;
A = 2.036646194262089E+08;
W = 1.619362528753911E+02*pi/180;
IN= 4.987675764908952E-01*pi/180;
TA= 1.932420048055015E+02*pi/180;

kep1 = [A EC IN OM W TA];
car1 = kep2car(kep1,astroConstants(4));

EC2= 1.780730530847942E-01;
OM2= 5.022651911496833E+01*pi/180;
A2 = 2.580534467305092E+08;
W2 = 2.444584674896587E+02*pi/180;
IN2= 3.836991178476995E+00*pi/180;
TA2= 2.962149957608689E+02*pi/180;

kep2 = [A2 EC2 IN2 OM2 W2 TA2];
car2 = kep2car(kep2,astroConstants(4));

kepm1=uplanet(6652,4)
carm1 = kep2car(kepm1,astroConstants(4));

kepm2=uplanet(6653,4)
carm2 = kep2car(kepm2,astroConstants(4));

figure()
hold on
grid on
plot3(car1(1),car1(2),car1(3),'x')
plot3(car2(1),car2(2),car2(3),'x')
plot3(kepm1(1),kepm1(2),kepm1(3),'o')
plot3(kepm2(1),kepm2(2),kepm2(3),'o')

%%
mu = astroConstants(4);
kep_E = uplanet(6653,3);
kep_M = uplanet(6653,4);
S_E = kep2car(kep_E,mu_s);  RE = S_E(1:3);
S_M = kep2car(kep_M,mu_s);  RM = S_M(1:3);

Rem = RM-RE;

figure()
hold on
grid on
