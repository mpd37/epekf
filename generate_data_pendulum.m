% Script to generate some data from a pendulum
%
% Copyright (C) 2012-2013 by 
% Marc Deisenroth and Shakir Mohamed
%
% Last modified: 2013-07-18

d  = 2;
e = 3;

FF =@(x,u) dynamics_pendulum_discretized(x,u);
%

sensor1 = [-2 0];
sensor2 = [2 0];
sensor3 = [-0.5 -0.5];

GG = @(x) [atan2( sin(x(2,:))-sensor1(2), cos(x(2,:))-sensor1(1) );
  atan2( sin(x(2,:))-sensor3(2), cos(x(2,:))-sensor3(1));
  x(1,:)
  ];



Q = diag([0.1, 0.3].^2);         % system noise
R = diag([0.1 0.05 0.1].^2);                 % measurement noise

m0 = zeros(d,1);
S0 = diag([0.5, pi/16].^2);

% generate time series
x(:,1) = gaussian(m0, S0);
y(:,1) = GG(x(:,1)) + chol(R)'*randn(e,1);
for t = 1:n-1
  u(:,t) = 4*rand(1)-2;
  x(:,t+1) = FF(x(:,t), u(:,t)) + chol(Q)'*randn(d,1);
  y(:,t+1) = GG(x(:,t+1)) + chol(R)'*randn(e,1);
end

u(:,t+1) = 0;


