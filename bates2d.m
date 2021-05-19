function main
    clc;
    clear all;
    close all;

    % Independent Variables
    r0 = 6697; 
    v0 = 10.8462;
    flightPathAngle0 = 0;
    arrivalAngle = 30;

    muEarth = 3.986012e5;  % Earth gravitational parameter
    muMoon = 4902.8; % Moon gravitational parameter 
    D = 384400; % radius of lunar orbit (km)
    Rs = 66300; % 66182;  % lunar sphere of influence (km)
    rm = 1738; % The radius of the Moon (km)

    E = v0^2/2 - muEarth/r0;
    h = r0*v0*cosd(flightPathAngle0);

    r1 = sqrt(D^2+Rs^2-2*D*Rs*cosd(arrivalAngle));
    if 2*(E+muEarth/r1) < 0
        fprintf("Insufficient Energy\n");
        return;
    else
        v1 = sqrt(2*(E+muEarth/r1));
    end
    
    flightPathAngle1 = acosd(h/r1*v1);
    gamma1 = asind((Rs/r1)*sind(arrivalAngle)); % degrees

% Determine geocentric trajectory elements to get the time of flight 
p = h^2/muEarth;
a = -muEarth/(2*E);
e = sqrt(1-p/a);
trueAnomaly0 = 0; % burn out at perigee, degrees
trueAnomaly1 = acosd((p-r1)/(r1*e)); % degrees
E0 = 0;
E1 = acosd((e+cosd(trueAnomaly1))/(1+e*cosd(trueAnomaly1))); % in degrees
timeOfFlight = sqrt(a^3/muEarth)*((deg2rad(E1)-e*sin(deg2rad(E1)))-(deg2rad(E0)-e*sin(deg2rad(E0))));
timeOfFlightHours = timeOfFlight/3600;
omegaMoon = 2.6847e-6; % rad/s
gamma0 = rad2deg(deg2rad(trueAnomaly1 - trueAnomaly0 - gamma1) - omegaMoon*timeOfFlight); % phase angle at departure in degrees

% At the Moon's sphere of influence 
vm = 1.018;


v2 = sqrt(v1^2+vm^2-(2*v1*vm*cosd(flightPathAngle1-gamma1)));

epsilon2 = asind(vm/v2*cosd(arrivalAngle) - v1/v2*cosd(arrivalAngle + gamma1 - flightPathAngle1));

% something not in Bates. TLI to Patch point yields Moon lead angle 
leadAngle = timeOfFlight * rad2deg(omegaMoon);

hMoon = Rs*v2*sind(epsilon2);

if hMoon <= 0
   retrograde = true;
else
   retrograde = false;
end

hMoon = abs(hMoon);

Em = v2^2/2 - muMoon/Rs;
pMoon = hMoon^2/muMoon;
eMoon = sqrt(1+ (2*Em*hMoon^2)/muMoon^2);

rp = pMoon/(1+eMoon);
hp = rp - rm;

vp = sqrt(2*(Em+(muMoon/rp)));


% fprintf("Initial Conditions\n");
% fprintf("------------------\n");
% fprintf("r0 = %0.1f km\n",r0);
% fprintf("v0 = %0.4f m/s\n", v0);
% fprintf("Initial Flight Path Angle = %0.0f deg\n", flightPathAngle0);
% fprintf("Lunar Arrival Angle = %0.0f deg\n", arrivalAngle);
% fprintf("\n");
% fprintf("Trajectory Conditions\n");
% fprintf("------------------\n");
% 
% fprintf("E = %0.4f km^2/s^2\n", E);
% fprintf("h = %0.0f km^2/s\n", h);
% fprintf("r1 = %0.0f km\n", r1);
% fprintf("v1 = %0.4f km/s\n", v1);
% fprintf("p = %0.0f km\n", p);
% fprintf("a = %0.0f km\n", a);
% fprintf("e = %0.4f \n", e);
% fprintf("Departure Phase Angle = %0.2f deg\n", gamma0);
% fprintf("\n");
% fprintf("Patch Conditions\n");
% fprintf("------------------\n");
% fprintf("Flight Path Angle = %0.3f deg\n", flightPathAngle1);
% fprintf("Phase Angle = %0.3f deg\n", gamma1);
% fprintf("True Anomaly = %0.2f deg\n", trueAnomaly1);
% fprintf("Eccentric Anomaly = %0.2f deg\n", E1);
% fprintf("Time of flight = %0.2f hours\n", timeOfFlightHours);
% fprintf("\n");
% fprintf("Moon Conditions\n");
% fprintf("------------------\n");
% fprintf("v2 = %0.4f km/s\n", v2);
% fprintf("Epsilon (Deviation) = %0.4f deg\n",epsilon2);
% fprintf("Eccentricity = %0.4f\n", eMoon);
% fprintf("h = %0.0f km^2/s\n", hMoon);
% fprintf("E = %0.4f \n", Em);
% if retrograde
%     fprintf("Retrograde\n");
% else
%     fprintf("Prograde\n");
% end
% fprintf("Lead angle = %0.2f deg\n", leadAngle);
% fprintf("Height at perilune: %0.1f km\n", hp);
% fprintf("Speed of probe relative to moon: %0.4f km/s\n", vp);
end