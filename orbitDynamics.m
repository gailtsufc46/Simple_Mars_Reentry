function udot = orbitDynamics(t,u,params,Atm)
%% Extract Parameters
mu = params.Mars.mu;
B = params.SC.B;
R = params.Mars.R;
%% Calculate Important State Characteristics
r = norm(u(1:3));
v = norm(u(4:6));
h = r-R;
[rho, ~, ~] = Atm(h);
global exit
if h<=5
    exit = 1;
end
%% Calculate State Derivatives
udot = zeros(length(u),1);
udot(1:3) = u(4:6);
udot(4:6) = -(mu*u(1:3)/r^3)-(1/2/B)*rho*v*u(4:6)*1000;
