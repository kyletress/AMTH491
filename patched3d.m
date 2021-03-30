% Kyle Tress
% AMTH 491 Senior Project
% Adapted from Orbital Mechanics for Engineering Students 

% Physical constants
Rs = 66183; % moon sphere of influence
Rm = 1737; % radius of the Moon in km
muEarth = 398600; % Earth gravitational parameter
muMoon = 4902.8; % Moon gravitational parameter 

% Independent Variables 
t = juliandate(datetime('2020-05-04 12:00:00'));
lambda = 50; % arrival angle 
r0 = 6698; % radius 
rightAscension = 40; % at TLI
declination = 10; % at TLI
flightPathAngle = 10; % at TLI

% Moon geocentric equatorial state vector from ephemeris
[rm,vm] = planetEphemeris(t, 'Earth','Moon');

% s, unit vector along the earth-moon radial 
s = rm/norm(rm);

% 4. The instantaneous angular velocity of the moon at SOI intercept
wm = cross(rm, vm)/norm(rm)^2;
angularVelocityMoon = norm(wm);

% 5. Geocentric position vector at TLI
r0 = [r0*(cosd(rightAscension)*cosd(declination)),r0*(sind(rightAscension)*cosd(declination)),r0*sind(declination)];

% 6. Unit normal to the plan of s and w1
w1 = cross(r0,rm)/norm(cross(r0,rm));

% 7. Calculate b, the unit normal to the plan of s and w1. (B-plane
% targeting in STK? Review. 
b = cross(w1,s)/norm(cross(w1,s));

% 8. Calculate n, the unit vector from the center of the moon to the patch
% point
n = (-cosd(lambda)*s)+sind(lambda)*b;

% 9. r2, the position vector of the patch point relative to moon 
r2 = Rs*n;

% 10. Position vector of the patch point relative to the Earth
r1 = rm + r2;

% 11. Sweep angle 
sweepAngle = acosd(dot(r0/norm(r0), r1/norm(r1)));

% 12. Angular momentum of the translunar trajectory 
h1 = sqrt(muEarth*norm(r0))*sqrt((1-cosd(sweepAngle))/(norm(r0)/norm(r1)+sind(sweepAngle)*tand(flightPathAngle)-cosd(sweepAngle)));

% 13. Lagrange Coefficients 
f = 1-((muEarth*norm(r1))/h1^2)*(1-cosd(sweepAngle));
g = ((norm(r0)*norm(r1))/h1)*sind(sweepAngle);
gdot = 1-((muEarth*norm(r0))/h1^2)*(1-cosd(sweepAngle));

% 14. Velocity at TLI and at the patch point 
v0 = 1/g*(r1-f*r0);
v1 = 1/g*(gdot*r1-r0);

% 15. Radial component of velocity at TLI
vr0 = dot(v0, r0/norm(r0));

% 16. Eccentricity vector - the norm must be less than 1!
e1 = (cross(v0, cross(r0,v0)))/muEarth - r0/norm(r0);
% norm(e1)

% 17. Semimajor axis and period of the translunar trajectory 
a1 = (norm(h1)^2/muEarth)*1/(1-norm(e1)^2);
T1 = 2*pi*sqrt(norm(a1)^3/muEarth);

% 18. Perifocal unit vectors p, q, and w (w was calculated above in step 6)
p1 = e1/norm(e1);
q1 = cross(w1,p1);

% 19. True anomoly at TLI. Note the velocity in step 16. Is it gte or lte 0
trueAnomoly = acosd(dot(p1,r0/norm(r0))); % REMINDER: this is in degrees, you'll need it in rads later later. 

% 20. Time since perigee at TLI 
E = 2*atan(sqrt((1-norm(e1))/(1+norm(e1)))*tan(deg2rad(trueAnomoly)/2));
meanAnomoly = E - norm(e1)*sin(E); % in radians
t0 = (meanAnomoly/(2*pi) * T1); % in seconds 

% 21. true anomoly at patch point 
theta1 = trueAnomoly + sweepAngle; % in degrees. 

% 22. Time since perigee at the patch point 
E1 = 2*atan(sqrt((1-norm(e1))/(1+norm(e1)))*tan(deg2rad(theta1)/2));
meanAnomoly1 = E1 - norm(e1)*sin(E1);
t1 = (meanAnomoly1/(2*pi) * T1); % in seconds 

% 23. flight time from TLI to patch point
flightTime = (t1 - t0)/3600; % in hours 


% II. Determine the lunar approach trajectory inside SOI (relative to the
% moon)

% 1. v2, velocity of spacecraft relative to the moon at the patch point 
v2 = v1 - vm;

% 2. Radial speed relative to the moon, at the patch point 
vr2 = dot(v2, r2/norm(r2));

% 3. Angular momentum of the trajectory 
h2 = cross(r2,v2);

% 4. Eccentricity vector of the trajectory relative to the moon. Must be
% greater than 1
e2 = (1/muMoon)*cross(v2, h2)-(r2/norm(r2));
% norm(e2)

% 5. Perilune radius and altitude of hyperbolic lunar approach 
rp2 = (norm(h2)^2/muMoon)*(1/(1+norm(e2)));
zp2 = rp2 - Rm;

% 6. Perifocal unit vectors of the hyperbolic approach trajectory 
p2 = e2/norm(e2);
w2 = h2/norm(h2);
q2 = cross(w2,p2);

% 7. orthogonal unit vectors directed along rotating xyz moon-fixed frame 
i = rm/norm(rm);
k = wm/norm(wm);
j = cross(k,i);

% 8. Instantaneous direction cosine matrix (transform from geocentric XYZ
% to moon fixed xyz
Q = [i(1) i(2) i(3); j(1) j(2) j(3); k(1) k(2) k(3)];

% 9. Transform components of the perifocal unit vectors into the moon-fixed xyz frame 
p2m = Q*p2';
q2m = Q*q2';

% 10. Time at the patch point
ur2 = r2/norm(r2);
theta2 = 360 - acosd(dot(p2,ur2));

F = 2*atanh(sqrt((norm(e2)-1)/(norm(e2)+1))*tan(deg2rad(theta2)/2));
meanAnomolyHyperbola = norm(e2)*sinh(F)-F;
t2 = meanAnomolyHyperbola * norm(h2)^3/(muMoon^2*(norm(e2)^2-1)^(3/2));




