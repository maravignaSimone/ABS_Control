clc
clear 
close all

%% DEFINE VARIABLES
%Params
syms Cd m R J real
g = 9.81;
%State 
syms v omega theta2 real
%Controls
syms Tb real
%Disturbances
syms theta1 theta3 slope wind dist_theta2 real
%Noise
syms nu_w nu_v real
%Reference
syms lambda_star real

%% DEFINE VECTORS
x = [v; omega; theta2];
u = Tb;
d = [theta1; theta3; wind; slope; dist_theta2];
nu = [nu_w; nu_v];
ref = lambda_star;
w = [d; nu; ref];

%% DEFINE EQUATIONS
% Aerodynamic drag
D = Cd*(v-wind)^2;

% Wheel load
N = (m*g)*cos(slope);

% Lambda (slip ratio when braking)
lambda = ((omega*R)-v)/v;

% Friction coefficient (Burckhardt model) (sign(lambda) = - 1 when braking)
mu = sign(lambda) * theta1 * (1-exp(-abs(lambda)*theta2))-(lambda*theta3);

% Friction forces
Fx = N * mu;

% System dynamics
f = [(Fx/m)-(D/m)-g*sin(slope); R*(-Fx/J)+Tb/J; dist_theta2];

% System outputs
h = [omega * (1 + nu_w); v + nu_v];
%he = (((omega*(1 + nu_w)*R)-(v + nu_v))/(omega*(1 + nu_w)*R)) - lambda_star;
he = (((omega*(1 + nu_w)*R)-(v + nu_v))/(v + nu_v)) - lambda_star;

%% computing matrices
A = jacobian(f, x);
B1 = jacobian(f, u);
B2 = jacobian(f, w);

C = jacobian(h, x);
D1 = jacobian(h, u);
D2 = jacobian(h, w);

Ce = jacobian(he, x);
D1e = jacobian(he, u);
D2e = jacobian(he, w);

%% defining matrices functions
matlabFunction(A, 'File', 'ABS_Amatrix');
matlabFunction(B1, 'File', 'ABS_B1matrix');
matlabFunction(B2, 'File', 'ABS_B2matrix');

matlabFunction(C, 'File', 'ABS_Cmatrix');
matlabFunction(D1, 'File', 'ABS_D1matrix');
matlabFunction(D2, 'File', 'ABS_D2matrix');

matlabFunction(Ce, 'File', 'ABS_Cematrix');
matlabFunction(D1e, 'File', 'ABS_D1ematrix');
matlabFunction(D2e, 'File', 'ABS_D2ematrix');
