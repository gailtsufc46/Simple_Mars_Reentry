%% Mars Re-entry Simulation
% Garrett Ailts
clear vars, close all

%% User Defined Parameters
h0 = 200; % Starting altitude in km
nOrbits = 2; % Number of orbits to simulate
delta_v = 25e-3; % Initial speed change from orbit velocity in km/s
B = 50; % Ballistic Coeffiecient of spacecraft in kg/m^2
d = 10; % Spacecraft diameter in m
%% Define Mars and Spacecraft Parameters
parameters.Mars.mu = 4.282837e4; % Gravitational parameter of Mars in km^3/s^2
parameters.Mars.R = 3396.2; % Equatorial radius of Mars in km
parameters.Mars.gamma = 1.68; % Adiabatic constant of the martian atmosphere, dimensionless
parameters.Mars.mu0 = 12.17e-6; % Dynamic viscosity of Mars in N-s/m^2
parameters.Mars.S = 222.2; % Effective temperature (Sutherlands constant) of the Martian atmosphere in degrees K
parameters.Mars.Mco2 = 0.044; % Molecular weight of CO2 atmosphere in kg/mol
parameters.Mars.Rgas = 8.134; % Universal gas constant in J/mol-K
parameters.SC.B = B; % Spacecraft ballistic coeffiecent in kg/m^2
rad2deg = 180/pi;
%% Create Necessary Function Handles
Atm = @(h) Mars_Atm(h,parameters);
Mach = @(v,s) v./s;
Re = @(v,nu,D) v*D./nu;
CF = @(M,Re) (0.65+0.399*((2/pi)*atan(10-M)+1))./sqrt(Re);
qavg = @(rho,v,CF) 0.25*rho.*v.^3.*CF;
%% Create Initial State Variables
v_parking = sqrt(parameters.Mars.mu/(parameters.Mars.R+h0));
T = 2*pi*sqrt((parameters.Mars.R+h0)^3/parameters.Mars.mu);
simTime = T*nOrbits;
t = linspace(0,simTime,simTime);
x0 = (h0 + parameters.Mars.R); y0 = 0; z0 = 0; vx0 = 0; vy0 = v_parking - delta_v; 
vz0 = 0; 
% Starting position in ECI frame.
u0 = [x0; % Initial x position in m
      y0; % Initial y position in m
      z0; % Initial z position in m
      vx0; % Intial x vecolicty in m/s
      vy0; % Intial y vecolicty in m/s
      vz0]; % Intial z vecolicty in m/s 
%% Run Orbit Simulation
uout = odeRK4(@(t,u) orbitDynamics(t,u,parameters,Atm),t,u0);
t = t(1:length(uout));
%% Preallocate Variables For Heating and Flight Path Angle Determination
rho = zeros(1,length(t));
h = zeros(1,length(t));
nu = zeros(1,length(t));
s = zeros(1,length(t));
r_mag = zeros(1,length(t));
v_mag = zeros(1,length(t));
fpa = zeros(1,length(t));
%% Obtain Flight Path Angle and Heating Rate At every Time Step
for i=1:length(t)
    r_mag = norm(uout(1:3,i));
    v_mag(i) = norm(uout(4:6,i));
    h(i)=r_mag-parameters.Mars.R; % Get altitude in km at every timestep
    ehat_up = uout(1:3,i)/r_mag;
    ehat_vel = uout(4:6,i)/v_mag(i);
    dp = dot(ehat_up,ehat_vel);
    alpha = acos(dp);
    fpa(i) = pi/2-alpha;
    [rho(i), nu(i), s(i)] = Atm(h(i));
end
Cf = CF(Mach(v_mag*1000,s),Re(v_mag*1000,nu,d));
q = qavg(rho,v_mag*1000,Cf);
%% Display Reentry Analysis
plot(t/3600,h);
xlabel('Time (hr)');
ylabel('Altitude (km)');
title('Orbit Altitude vs Time');
[qmax, qidx] = max(q);
h_qmax = h(qidx);
v_final = v_mag(end);
fpa_final = fpa(end)*rad2deg;
duration = t(end);
heating_max = sprintf('The max heating rate is %f W and occurs at an alitude of %f km',qmax,h_qmax);
trajectory_final = sprintf('The final velocity is %f km/s at a flight path angle of %f degrees',v_final,fpa_final);
sim_duration = sprintf('The re-entry time was %f s',duration);
disp(heating_max);
disp(trajectory_final);
disp(sim_duration);
