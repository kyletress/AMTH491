% Kyle Tress
% Applied Math 491, Senior Project 
% This script analyzes 72 million patched conic lunar trajectories and plots 
% those with a perilune altitude between 499.5 and 500km. 

function main
    clc; clear all; close all;
    
    % Constants
    muEarth = 3.986012e5;  % Earth gravitational parameter
    muMoon = 4902.801; % Moon gravitational parameter 
    D = 384400; % radius of lunar orbit (km)
    Rs = 66300; % 66182;  % lunar sphere of influence (km)
    rm = 1737.1; % The radius of the Moon (km)
    vm = 1.018; % moon's orbital velocity 
    
    % Initial position and orbital velocity
    altitude = 300;
    r1 = altitude + 6377;
    v0 = sqrt(muEarth/r1);
    deltaV = 3.131;
    thrustAngle = 0;
    arrivalAngle = 46;

     v1 = sqrt(v0^2+deltaV^2+(2*v0*deltaV*cosd(thrustAngle)));
     gamma1 = acosd((deltaV^2-v0^2-v1^2)/(-2*v0*v1))*sign(thrustAngle);
        E = v1^2/2 - muEarth/r1;
        h = r1*v1*cosd(gamma1);
        r2 = sqrt(D^2+Rs^2-(2*D*Rs*cosd(arrivalAngle)));

        if 2*(E+muEarth/r2) < 0
            fprintf("Insufficient energy");
            return; % insufficient energy to reach the Moon
        else
            v2 = sqrt(2*(E+muEarth/r2));
        end

        gamma2 = acosd(h/(r2*v2));
        beta2 = asind(Rs/r2*sind(arrivalAngle));
        v3 = sqrt(v2^2+vm^2-2*v2*vm*cosd(gamma2-beta2));
        epsilon3 = asind(vm/v3*cosd(arrivalAngle) - v2/v3*cosd(arrivalAngle + beta2 - gamma2));

        % At the patch point
        hMoon = Rs*v3*sind(epsilon3);
        Em = v3^2/2 - muMoon/Rs;
        pMoon = hMoon^2/muMoon;
        eMoon = sqrt(1+ (2*Em*hMoon^2)/muMoon^2);
        rp = pMoon/(1+eMoon);
        hp = rp - rm;
        vp = sqrt(2*(Em+(muMoon/rp)));

        p = h^2/muEarth;
        a = -muEarth/(2*E);
        e = sqrt(1-p/a);
        trueAnomaly1 = acosd((p-r1)/(r1*e));
        trueAnomaly2 = acosd((p-r2)/(r2*e));
        E1 = acosd((e+cosd(trueAnomaly1))/(1+e*cosd(trueAnomaly1))); 
        E2 = acosd((e+cosd(trueAnomaly2))/(1+e*cosd(trueAnomaly2))); 
        timeOfFlight = sqrt(a^3/muEarth)*((deg2rad(E2)-e*sin(deg2rad(E2)))-(deg2rad(E1)-e*sin(deg2rad(E1))));
        omegaMoon = 2.6847e-6; % rad/s
        omega1 = trueAnomaly2 - trueAnomaly1 - beta2 - rad2deg(omegaMoon*timeOfFlight);
        fprintf("%0.4f  %0.4f  %0.0f %0.2f  %0.1f  %0.1f  %0.1f %0.1f %0.1f\n",deltaV, epsilon3, a, omega1, arrivalAngle, thrustAngle, timeOfFlight/3600, hp, trueAnomaly2);
end