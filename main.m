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
ld = 5; % disturbance
r = ld+q+lm; % exogenous

%% nominal parameters
%LAMBORGHINI HURACAN
% params
m = 387.5; % [kg] quarter-vehicle mass
R_w = 0.33; % [m] wheel radius
Ji = 0.45; % [kg*m^2] wheel inertia
% aerodynamics
rho = 1.225; % [kg/m^3] air density
S = 1.75; % [m^2] cross section
CD = 0.39; % [-] drag resistance coefficient
Cd = 0.5*rho*S*CD; % drag coefficient (1/2)*rho*S*Cd

%Other vehicle
% % params
% m = 250; % [kg] quarter-vehicle mass
% R_w = 0.3; % [m] wheel radius
% Ji = 1; % [kg*m^2] wheel inertia
% % aerodynamics
% rho = 1.225; % [kg/m^3] air density
% S = 1.6*1.2; % [m^2] cross section
% CD = 2.2; % [-] drag resistance coefficient
% Cd = 0.5*rho*S*CD; % drag coefficient (1/2)*rho*S*Cd

%% Linearization Conditions
slope0 = 0;
wind0 = 0;
nu_w = 0;
nu_v = 0;
theta1=1.28;
theta3=0.52;
dist_theta2_0 = 0;

% friction coefficients
%dry asphalt
theta0_2=23.99;
lambda_star = -0.17;
%wet asphalt
%theta1 = 0.86;
%theta0_2 = 33.82;
%theta3 = 0.35;
%lambda_star= -0.131;

v0 = 130/3.6; % [m/s] initial velocity

%%

% function for finding lambda0 of the equilibrium point
eqn1 = g*cos(slope0)*(sign(lambda_eq) * theta1 * (1-exp(-abs(lambda_eq)*theta0_2))-(lambda_eq*theta3)) - (Cd/m)*(v0-wind0)^2 -g*sin(slope0);
eqn2 = Tb_eq/Ji - (R_w/Ji)*(m*g*cos(slope0)*(sign(lambda_eq) * theta1 * (1-exp(-abs(lambda_eq)*theta0_2))-(lambda_eq*theta3)));

f1 = matlabFunction(eqn1, 'Vars', [Tb_eq, lambda_eq]);
f2 = matlabFunction(eqn2, 'Vars', [Tb_eq, lambda_eq]);

fun = @(var) [f1(var(1), var(2)); f2(var(1), var(2))];

var0 = [0; 0];

var = fsolve(fun, var0);

Tb_eq = var(1);
lambda_eq = var(2);

%omega0 = (lambda_eq*v0 + v0)/R_w; %braking
omega0 = - v0/(R_w*(lambda_eq-1)); %driving

%% Linearized plant
%state
Amat = ABS_Amatrix(Cd,Ji,R_w,m,omega0,slope0,theta1,theta0_2,theta3,v0,wind0);
B1mat = ABS_B1matrix(Ji);
B2mat = ABS_B2matrix(Cd,Ji,R_w,m,omega0,slope0,theta1,theta0_2,theta3,v0,wind0);

%output
Cmat = ABS_Cmatrix(nu_w);
D1mat = ABS_D1matrix;
D2mat = ABS_D2matrix(omega0);

%error
% he = omega - omega_ref
% Cemat = ABS_Cematrix(R_w,lambda_star,nu_w);
% D1emat = ABS_D1ematrix;
% D2emat = ABS_D2ematrix(R_w,lambda_star,omega0);

% he: lambda - lambda_star
Cemat = ABS_Cematrix(R_w,nu_v,nu_w, omega0, v0);
D1emat = ABS_D1ematrix;
D2emat = ABS_D2ematrix(R_w,nu_v,nu_w,omega0,v0);
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
tb0 = Tb_eq;

x0 = [v0; omega0; theta0_2];

u0 = tb0;

d0 = [theta1; theta3; wind0; slope0; dist_theta2_0];

w0 = [d0; nu_w; nu_v; lambda_star];

y0 = [x0(2); x0(1)];

e0 = 0;

v_init = v0; % [m/s] vehicle speed
omega_init = v_init/R_w; % [rad/s] rear wheel speed

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

eps1max = v0;
eps2max = omega0;
eps3max = 0.0005;

Q = inv(length(Ceps)*diag([eps1max ^2,eps2max^2, eps3max^2]));

umax = 700;

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

%% OBSERVER
lambda_d = 0;
Ad = Amat.';
Bd = Cmat.';
Cd = B2mat.';
Dd = D2mat.';

q1 = 0.05;
q2 = 0.05;
q3 = 10;
q4 = atan(0.1);
q5 = 0.1;
q6 = 0;
q7 = 0;
q8 = 0;

Qd = diag([q1^2, q2^2, q3^2, q4^4, q5^2, q6^2, q7^2, q8^2]);

std_tacho = 0.05;
std_GPS = 0.01;

Rd = diag([std_tacho^2, std_GPS^2]);

barRd = Rd + Dd.'*Qd*Dd;

Amd = Ad + lambda_d*eye(n);
Bmd = Bd;
Qmd = Cd.'*Qd*Cd;
Rmd = barRd;
Smd = Cd.'*Qd*Dd;
Emd = eye(n);
Gmd = 0;

[X, Kmd, L] = icare (Amd, Bmd, Qmd, Rmd, Smd, Emd, Gmd);
KO = Kmd.';

AO = Amat - KO*Cmat;
BO = B1mat - KO*D1mat;
CO = eye(n);
DO = zeros(n, q);
XOinit = x_init - x0;

%% RUN THE SIMULATOR
PLANT = 1; % 0 = linear, 1 = nonlinear
TimeSpan = 6;
DT = 1e-4;

% plots
out = sim('SimulinkModel',TimeSpan);
f = figure('Name','Measured Data','NumberTitle','off');
t=tiledlayout(2, 2);
nexttile;
plot(out.linear.y, LineWidth=1);
legend("omega", "v");
title("Measurements $\tilde{y}$", Interpreter="latex");
nexttile;
plot(out.linear.e, LineWidth=1);
title("Error $\tilde{e}$", Interpreter="latex");
nexttile;
plot(out.linear.u, LineWidth=1);
title("Control: Torque $\tilde{u}$", Interpreter="latex");
nexttile;
plot(out.lambda, LineWidth=1);
title("Slip Ratio $\lambda$", Interpreter="latex");
title(t, strcat("Result of the control with $\epsilon_1$=",num2str(eps1max), ", $\epsilon_2$=", num2str(eps2max), ", $\epsilon_3$=", num2str(eps3max), ", $u_{max}$=", num2str(umax)), Interpreter="latex");