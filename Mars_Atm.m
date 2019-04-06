function [rho, nu, s] = Mars_Atm(h,params)
% Function that returns density, temperature, and kinematic viscosity at a 
% specific altitude in the Martian atmosphere. 
%
% Returns the density in kg/m^3, the temperature in degrees C, and the
% kinemtaic viscosity in m^2/s
%
%% Define Parameters and Extract parameters
h_m = h*1000;
a0 = 49.8118119899434;
a1 = -5.9123700325916;
a2 = -3.5638800977374;
a3 = 0.380908561109888;
gamma = params.Mars.gamma;
mu0 = params.Mars.mu0;
Rgas = params.Mars.Rgas;
S = params.Mars.S;
Mco2 = params.Mars.Mco2;
%% Define levels of the Martian Atmosphere
if h<7
    T=-31-0.000998.*h_m;
elseif h>=7 && h<=65
    T = -23.4-0.00222.*h_m;
end
P = 0.699.*exp(-0.00009.*h_m);
if h>65
   T = -167.7;
   rho = 0.88325.*exp(a0+a1*log(h)+a2*log(h).^2+a3*log(h).^3);
else
   rho = P./(.1921.*(T+273.1)); % Density in kg/m^3 
end
T = T+273.15;
mu = mu0*(T/242.15)^(3/2)*((242.15+S)/(T+S));
nu = mu/rho; % Kinematic viscosity in m^2/s;
s = sqrt(gamma*Rgas*T/Mco2); % Speed of sound in m/s

