function [dot_x,y,e] = F(x,u,w)

%% VARIABLES

% STATES
v = x(1); % [m/s] vehicle speed
omega = x(2); % [rad/s] wheel speed
theta2 = x(3); % theta2

% DISTURBANCE
wind = w(1);  % [m/s] wind speed
slope = w(2); % [rad] road slope

% NOISE
nu = w(3:4); % [GNSS sensor [m/s]
%               tone-wheel [rad/s]];

% REFERENCE
ref = w(5);

% CONTROL
Tb = u; % [Nm] rear wheel torque

%% ACTUAL PLANT PARAMETERS

m = 250; % [kg] quarter-vehicle mass
g = 9.81; % [m/s^2] gravity acceleration
Ji = 1; % [kg*m^2] rear wheel inertia
R = 0.3; % [m] rear wheel radius
theta1=1.28;
theta3=0.52; 
%% STATE DYNAMICS

% slip ratio
lambda_w = lambda(v,omega*R);
theta = [theta1 theta2 theta3];
mu_w = mu(lambda_w,theta);

% wheel load
N = m*g*cos(slope);

% friction forces
Fx = N * mu_w;

% system dynamics
dot_x = [-g*sin(slope)+Fx/m-D(v-wind)/m; (Tb - Fx*R)/Ji; 0];

%% MEASUREMENTS

y = [omega*(1+nu(1)); v + nu(2)];

%% CONTROLLED VARIABLES

lambda_star = ref;
e = y(1)-(y(2)/R)*(1+lambda_star);
% e = ((y(1)*R)-y(2))/y(2) - lambda_star;
end

%% SUPPORT FUNCTIONS

% slip ratio
function L = lambda(v,omegar)
if abs(v-omegar) <= 1e-4  % pure rolling / at rest
    L = 0;
else
    vmax = max(abs(v),abs(omegar));
    vmax = max(vmax,abs(omegar-v));
    L = (omegar-v)/vmax; % driving / braking
end
end

% friction model
function m = mu(a1,Theta)
m = sign(a1)*Theta(1)*(1-exp(-abs(a1)*Theta(2)))-a1*Theta(3); % Burckhardt model
end


% air drag
function out = D(v)
rho = 1.225; % [kg/m^3] air density
S = 1.6*1.2; % [m^2] cross section
CD = 2.2; % [-] drag resistance coefficeint

out = 0.5*rho*S*CD*(v^2);
end