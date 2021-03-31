% Kyle Tress
clc;
radiusEarth = 6378;
orbitalAltitude = 320;
alpha0 = 28;
arrivalAngle = 55;
flightPathAngle = 6;
Rs = 66183; % lunar sphere of influence 
muEarth = 398600; % Earth gravitational parameter
muMoon = 4902.8; % Moon gravitational parameter 
D = 384400; % assume circular orbit for the Moon 
velocityMoon = 1.0183;

% Initial position vector
r0 = radiusEarth + orbitalAltitude;
r0 = [-r0*cosd(alpha0) -r0*sind(alpha0) 0];
ur0 = r0/norm(r0);

r2 = Rs;
r2 = [-Rs*cosd(arrivalAngle) Rs*sind(arrivalAngle) 0];
ur2 = r2/norm(r2);

rm = [384400 0 0];
r1 = rm + r2;
ur1 = r1/norm(r1);

sweepAngle = acosd(dot(ur0,ur1));

h1 = sqrt(muEarth*norm(r0))* sqrt((1-cosd(sweepAngle))/(norm(r0)/norm(r1)+sind(sweepAngle)*tand(flightPathAngle)-cosd(sweepAngle)));

% Lagrange Coefficients 
f = 1-((muEarth*norm(r1))/h1^2)*(1-cosd(sweepAngle));
g = ((norm(r0)*norm(r1))/h1)*sind(sweepAngle);
gdot = 1-((muEarth*norm(r0))/h1^2)*(1-cosd(sweepAngle));

v0 = 1/g * (r1 - f*r0);
norm(v0);
vr0 = dot(v0,ur0); % The radial component of the velocity at v0

v1 = 1/g * (gdot*r1 - r0);
norm(v1);
vr1 = dot(v1,ur1);

e1 = (1/muEarth)*((norm(v0)^2-(muEarth/norm(r0)))*r0-norm(r0)*norm(vr0)*v0);
norm(e1);

% Perifocal unit vectors 
p1 = e1/norm(e1);
w1 = cross(r0,v0)/norm(h1);
q1 = cross(w1,p1);

% Semimajor axis
a1 = norm(h1)^2/muEarth * 1/(1-norm(e1)^2);
% Period
T1 = 2*pi * sqrt(a1^3/muEarth);
% True anomoly of the injection point
trueAnomolyTLI = acosd(dot(p1,ur0));

% Time at TLI 
E = 2*atan(sqrt((1-norm(e1))/(1+norm(e1)))*tan(deg2rad(trueAnomolyTLI)/2));
meanAnomoly = E - norm(e1)*sin(E); % in radians
t0 = (meanAnomoly/(2*pi) * T1); % in seconds 

trueAnomolyPatch = trueAnomolyTLI + sweepAngle;
E = 2*atan(sqrt((1-norm(e1))/(1+norm(e1)))*tan(deg2rad(trueAnomolyPatch)/2));
meanAnomoly = E - norm(e1)*sin(E); % in radians
t1 = (meanAnomoly/(2*pi) * T1);

flightTime = t1 - t0;

vm = [0 1.0183 0];
v2 = v1 - vm;
angularVelocityMoon = velocityMoon/D;

leadAngle = rad2deg(angularVelocityMoon)*flightTime;


fprintf(' Output\n\n');
fprintf(' Orbital Parameters\n');
fprintf(' --------------------------------------------\n');
fprintf(' Eccentricity:            %0.4f\n', norm(e1));
fprintf(' Semimajor Axis:          %0.0f km\n', a1);
fprintf(' Period:                  %0.2f days\n', T1/(3600*24));
fprintf(' True Anomoly at TLI:     %0.4f deg\n', trueAnomolyTLI);
fprintf(' Flight Time:             %0.4f hrs\n', flightTime/3600);
fprintf(' Lead Angle:              %0.2f deg\n', leadAngle);
fprintf('\n');
fprintf(' Spacecraft Parameters\n');
fprintf(' --------------------------------------------\n');
fprintf(' Velocity Relative to Moon:   %0.5f km/s\n', norm(v2));
