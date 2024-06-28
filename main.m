clc
clear
close all

symbolic;
syms Tb_eq lambda_eq real

%% DECLARATION OF VECTOR DIMENSIONS

n = 3; % state
p = 1; % control
q = 2; % measurement
lm = 1; % regulated output
ld = 2; % disturbance
r = ld+q+lm; % exogenous

%% nominal parameters

% params
m = 250; % [kg] quarter-vehicle mass
R = 0.3; % [m] wheel radius
Ji = 1; % [kg*m^2] wheel inertia
theta1=1.28;
theta3=0.52;
% aerodynamics
rho = 1.225; % [kg/m^3] air density
S = 1.6*1.2; % [m^2] cross section
CD = 2.2; % [-] drag resistance coefficient
Cd = 0.5*rho*S*CD; % drag coefficient (1/2)*rho*S*Cd

%% Linearization Conditions
slope0 = 0;
wind0 = 0;
nu_w = 0;
nu_v = 0;

% friction coefficients
%dry asphalt
theta0_2=23.99;
lambda_star = -0.17;
%wet asphalt
%theta1 = 0.86;
%theta0_2 = 33.82;
%theta3 = 0.35;
%lambda_star= -0.131;

v0 = 14; % [m/s] initial velocity

%%

% function for finding lambda0 of the equilibrium point
eqn1 = g*cos(slope0)*(sign(lambda_eq) * theta1 * (1-exp(-abs(lambda_eq)*theta0_2))-(lambda_eq*theta3)) - (Cd/m)*(v0-wind0)^2 -g*sin(slope0);
eqn2 = Tb_eq/Ji - (R/Ji)*(m*g*cos(slope0)*(sign(lambda_eq) * theta1 * (1-exp(-abs(lambda_eq)*theta0_2))-(lambda_eq*theta3)));

f1 = matlabFunction(eqn1, 'Vars', [Tb_eq, lambda_eq]);
f2 = matlabFunction(eqn2, 'Vars', [Tb_eq, lambda_eq]);

fun = @(var) [f1(var(1), var(2)); f2(var(1), var(2))];

var0 = [0; 0];

var = fsolve(fun, var0);

Tb_eq = var(1);
lambda_eq = var(2);

%omega0 = (lambda_eq*v0 + v0)/R;
omega0 = - v0/(R*(lambda_eq-1));

%% Linearized plant
%state
Amat = ABS_Amatrix(Cd,Ji,R,m,omega0,slope0,theta1,theta0_2,theta3,v0,wind0);
B1mat = ABS_B1matrix(Ji);
B2mat = ABS_B2matrix(Cd,Ji,R,m,omega0,slope0,theta1,theta0_2,theta3,v0,wind0);

%output
Cmat = ABS_Cmatrix(nu_w);
D1mat = ABS_D1matrix;
D2mat = ABS_D2matrix(omega0);

%error
Cemat = ABS_Cematrix(R,lambda_star,nu_w);
D1emat = ABS_D1ematrix;
D2emat = ABS_D2ematrix(R,lambda_star,nu_v,omega0,v0);

%% Jordan Canonical Form
[V,Vn,J] = JCF(Amat);

%% spider plot
n = length(Amat);
for kk = 1:n
    VV = Vn(kk,:);
    figure
    h = spider_plot(VV,...
        'AxesLimits', [zeros(1,n); 100*ones(1,n)],...
        'AxesLabels', {'$z_1$', '$z_2$','$z_3$'},...
        'AxesLabelsEdge', 'none',...
        'AxesInterpreter', 'latex',...
        'AxesTickLabels', 'data',...
        'Direction', 'counterclockwise',...
        'AxesDisplay', 'none',...
        'LabelFontSize', 24,...
        'AxesOffset', 0, ...
        'Color',0*ones(1,3), ...
        'Marker', {'o'},...
        'MarkerSize', 140, ...
        'MarkerFaceColor',{'w'},...
        'MarkerEdgeColor',{'k'});
    title(['$\tilde{x}_',num2str(kk),'$ modal components'],'Interpreter','latex','FontSize',24)
end

%% Reachability and Observability
Rmat = ctrb(Amat,B1mat);
Omat = obsv(Amat,Cmat);

disp(strcat("The reachable states are: ", num2str(rank(Rmat))));
disp(strcat("The observable states are: ", num2str(rank(Omat))));

%% INITIAL CONDITIONS
tb0 = -Tb_eq;

theta0 = [theta1; theta0_2; theta3];

x0 = [v0; omega0; theta0(2)];

u0 = tb0;

d0 = [wind0; slope0];

w0 = [d0; nu_w; nu_v; lambda_star];

y0 = [x0(2); x0(1)];

e0 = 0;

v_init = v0; % [m/s] vehicle speed
omega_init = v_init/R; % [rad/s] rear wheel speed

x_init = [v_init; omega_init; theta0_2];

%% KALMAN

[A_bar, B1_bar, C_bar, T, k] = ctrbf(Amat, B1mat, Cmat);

% Controllable part
Acc = A_bar(k+1:end, k+1:end);
Bcc = B1_bar(k+1:end, :);
Ccc = C_bar(:, k+1:end);

% Unontrollable part
Auc = A_bar(1:k, 1:k);
Buc = B1_bar(1:k, :);
Cuc = C_bar(:, 1:k);

% Augment with integral action
n_c = size(Acc, 1);

% Pick the equivalent state and tell what are the reachable ones
eq = T*x; disp(strcat("A reachable states is: ", string(eq(k+1:end))));

%% STATE FEEDBACK CONTROL + INTEGRAL ACTION

Ae = [Acc zeros(n_c, lm)
    Cemat(:, 1:k+1) zeros(lm, lm)];

Be = [Bcc;
    D1emat];

Ceps = eye(n_c+lm);

Deps = zeros(n_c+lm, p);

eps1max = 100;
eps2max = 100;
eps3max = 1;

Q = inv(length(Ceps)*diag([eps1max ^2,eps2max^2, eps3max^2]));

umax = 500;

R = inv(p*diag(umax ^2));

barR = R + Deps.'*Q*Deps;

alpha = 0;

Am = Ae + alpha*eye(n_c+lm);
Em = eye(n_c+lm);
Bm = Be;
Gm = 0;
Qm = Ceps.'*Q*Ceps;
Sm = Ceps.'*Q*Deps;
Rm = barR;

[X, Km, L] = icare (Am, Bm, Qm, Rm, Sm, Em, Gm);
K = -Km;
% Extract Ks and Ki from K
KS = K(:, 1:n_c);
KI = K(:, n_c+1:end);
%% RUN THE SIMULATOR
PLANT = 0; % 0 = linear, 1 = nonlinear
TimeSpan = 4;
DT = 1e-6;
%% Simulink model
%out = sim('SimulinkModel',TimeSpan);
%save CurrentWorkspace

%% PLOT RESULTS

%run plotResults