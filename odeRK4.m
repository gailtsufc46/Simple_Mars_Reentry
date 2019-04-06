function x  = odeRK4(rhs,time,x0)
% By Garrett Ailts
%
% MATLAB function that takes in a first order ODE representing the change
% of state of a dynamical system, a vector of discrete time, and an initial 
% state x0 as a column vector. The function then returns the state at every
% time step using a 4th order Runge-Kutta numerical integration method
%
% Function can except a state vector of any length, provided the passed 
% change of state equation can operate on all the state variables
%
%% Initialize State
x = zeros(length(x0),length(time));
x(:,1) = x0;
global exit;
exit = 0;
%% Propagate State Using Runge-Kutta Method
for i=1:length(time)-1
    h = time(i+1)-time(i);
    K1 = rhs(time(i),x(:,i));
    K2 = rhs(time(i)+h/2,x(:,i)+K1*h/2);
    K3 = rhs(time(i)+h/2,x(:,i)+K2*h/2);
    K4 = rhs(time(i+1),x(:,i)+K3*h);
    K = K1/6+K2/3+K3/3+K4/6;
    x(:,i+1) = x(:,i) + K*h;
    if exit==1
        x = x(:,1:i);
        break
    end
end
