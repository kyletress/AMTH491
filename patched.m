% Kyle Tress
% AMTH 491
% Patched Conic Approximation given R0, V0, lambda_1, and flight path angle

function main
% Earth-Moon properties 
mu = 398600; % Earth gravitational parameter
mu_moon = 4902.8;
Rs = 66300; % moon sphere of influence (km)
D = 384400; % km

% Initial Conditions 
r0 = 6697; % Injection altitude, km
v0 = 10.8462; % Velocity, km/s
flightPathAngle = 0;
lambda1 = 30; % arrival angle, degrees

% We first determine the energy and angular momentum of the orbit 

E = (v0^2/2) - (mu/r0); % energy
h = r0*v0*cosd(flightPathAngle); % angular momentum 

% Next we determine the radius at lunar arrival 
r1 = sqrt(D^2+Rs^2-2*D*Rs*cosd(lambda1));
v1 = sqrt(2*(E+mu/r1));
flightPathArrival = acosd(h/r1*v1);

gamma1 = asin(Rs/r1*sind(lambda1)); % Phase angle of the moon at arrival

% Next we need to calculate the line of flight to the moon 
% starting with the parameters p, a, e E0, E1

p = (h^2)/mu;
a = (-mu)/(2*E); % semimajor axis, km
e = sqrt(1-(p/a)); % eccentricity

% True anomolies 
w0 = acos((p-r0)/(r0*e)); 
w1 = acos((p-r1)/(r1*e));

% Eccentric anomolies 
E0 = acos((e+cos(w0))/(1+e*cos(w0)));
E1 = acos((e+cos(w1))/(1+e*cos(w1)));

% Now calculate the time of flight 
timeOfFlight = (sqrt((a^3)/mu)*((E1-(e*sin(E1))) - (E0-(e*sin(E0)))))/3600

% Phase angle at departure 
gamma0 = w1-w0-gamma1-((2.649e-6)*(timeOfFlight))

fprintf(" p: %2.0fkm\n",p);
fprintf(" a: %2.0fkm\n",a);
fprintf(" e: %2.4f\n",e);
fprintf(" Time of Flight: %2.4f hrs\n",timeOfFlight);
fprintf(" Phase Angle of Moon: %2.4f rads\n",gamma1);
fprintf(" True Anamolies: %2.3f, %2.3f\n", w0,w1);
fprintf(" Eccentric Anamolies: %2.3f, %2.3f\n", E0,E1);

t1 = datetime('2020-05-04 12:00:00');
d = juliandate(t1)

end 