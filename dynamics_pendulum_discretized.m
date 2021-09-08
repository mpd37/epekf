function dz = dynamics_pendulum_discretized(z,u)

% pendulum dynamics derived from Lagrangian
%
% state is given by [dtheta, theta]
%
% theta:   angular velocity
% dtheta:  angle
%
% (C) Copyright 2006-2013 
% by Marc Deisenroth, 2011-02-18

dt = 0.2;

m = 1;        % [kg]     mass of pendulum
l = 1;        % [m]      length of pendulum
g = 9.82;     % [m/s^2]  acceleration of gravity
I = m*l^2/12; % moment of inertia around midpoint of the pendulum

supersample = 2000;
dt2 = dt./supersample;

dz = z;

for t = 1:supersample
  f = (u-m*g*l*sin(dz(2,:))/2)./(I+m*l^2/4); % angular velocity
  dz = dz + dt2.*[f; dz(1,:)]; % anglular velocity/angle
end




