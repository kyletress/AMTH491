% Kyle Tress
% AMTH491 Senior Project 
% Simple Lunar Trajectories

function main

moon_apoapsis = 405500;
moon_periapsis = 363400;
semiMajorAxisMoon = semiMajorAxis(moon_periapsis,moon_apoapsis); % the semimajor axis of the moon
D = semiMajorAxisMoon; % Simplifying assumption, the moon is in a circular orbit around the earth with diameter D
mu = 398600; % Earth gravitational parameter 
v_m = orbitalSpeed(mu,semiMajorAxisMoon); % orbital speed of the moon
alt = 320; % the altitude of the spacecraft 
earthRadius = 6378; % radius of the Earth in km 
r_c = earthRadius + alt; 
v_c = orbitalSpeed(mu,r_c); % The orbital speed of the vehicle in parking orbit

% Hohmann Transfer Parameters
% Perigee is tangent to the circular orbit of the vehicle. 
% Apogee is tangent to the moon's orbit 
r_p = r_c;
r_a = semiMajorAxisMoon;
semiMajorAxisTransfer = semiMajorAxis(r_p,r_a);
eTransfer = eccentricity(r_a,r_p);
periodTransfer = period(semiMajorAxisTransfer) / 3600; % the period, converted to hours
timeOfFlight = periodTransfer/2 / 24; % in days 
h = angularMomentum(r_a,r_p); % the angular momentum of the transfer ellipse

% Now that we have the angular momentum, we can find the speeds at perigee
% and apogee (km/s)
v_p = h/r_p;
v_a = h/r_a;

% Delta-v requirements, km/s
deltav_p = v_p - v_c; % delta-v requiree to transfer from parking orbit to Hohmann transfer trajectory perigee
deltav_a = v_m - v_a; % delta-v required to transfer from Hohmann ellipse to moon's orbit
delta_v = deltav_p + deltav_a; % total delta-v requirement

fprintf(" Initial Parking Orbit Altitude: %2.2f km\n Time of Flight: %2.2f days\n Total Delta-v required: %2.4f km/s\n\n", alt,timeOfFlight, delta_v);

% Patched Conic Lunar Trajectories 
fprintf("Patched Conic Lunar Trajectories\n\n")
% Sphere of influence of the moon
Rs = 66183; % km 

r_m = [D,0,0]; % position vector of the moon's orbit, assumed to be constant magnitude 
v_m = [0,v_m,0];
alpha_0 = 28; % Angular position relative to the Earth-Moon line
gamma_0 = 6; % flight path angle 
lambda_0 = 55; % lunar arrival angle 
r_0 = [-r_c*cosd(alpha_0),-r_c*sind(alpha_0),0] % the position vector of the spacecraft at TLI 

% calculate the unit vector 
u_r0 = r_0/norm(r_0);

r_2 = [-Rs*cosd(lambda_0),Rs*sind(lambda_0),0]; % vector from the moon to the spacecraft at the patch point
% calculate the unit vector 
u_r2 = r_2/norm(r_2)

r_1 = r_m + r_2
norm(r_1) % magnitude of the distance to the patch point from the center of the earth 


    function answer = semiMajorAxis(a, p)
        % This function calculates the semimajor axis given the periapsis
        % and apoapsis
        answer =(a + p)/2;
    end

    function answer = orbitalSpeed(mu,D)
        % This function calculates the orbital speed given mu and the
        % diamter. Assumes circular orbit. km/s
        answer = sqrt(mu/D);
    end

    function answer = eccentricity(r_a, r_p)
        answer = (r_a-r_p)/(r_a+r_p);
    end

    function answer = period(a)
        % Returns the period in seconds
        answer = (2*pi*a^(3/2))/sqrt(mu);
    end

    function answer = angularMomentum(r_a,r_p)
        % Returns the angular momentum, km^2/s
        answer = sqrt(2*mu)*sqrt((r_a*r_p)/(r_a+r_p));
    end

end
