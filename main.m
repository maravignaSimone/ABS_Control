clc
clear
close all

symbolic;
syms lambda0 real

%% DECLARATION OF VECTOR DIMENSIONS

n = 5; % state
p = 1; % control
q = 2; % measurement
lm = 1; % regulated output
ld = 2; % disturbance
r = ld+q; % exogenous

%% nominal parameters

% params
m = 250; % [kg] quarter-vehicle mass
R = 0.3; % [m] wheel radius
Ji = 1; % [kg*m^2] wheel inertia

% aerodynamics
rho = 1.225; % [kg/m^3] air density
S = 1.6*1.2; % [m^2] cross section
CD = 2.2; % [-] drag resistance coefficient
Cd = 0.5*rho*S*CD; % drag coefficient (1/2)*rho*S*Cd

%% Linearization Conditions
v0 = 28; % [m/s] initial velocity

slope0 = 0;
wind0 = 0;
nu_w = 0;
nu_v = 0;

% friction coefficients
%dry asphalt
theta0_1=1.28;
theta0_2=23.99;
theta0_3=0.52;
lambda_star = -0.17;
%wet asphalt
%theta0_1 = 0.86;
%theta0_2 = 33.82;
%theta0_3 = 0.35;
%lambda_star= -0.131;

% function for finding lambda0 of the equilibrium point
l = (Cd*(v0-wind0)^2) + m*g*sign(lambda0)*theta0_1*(1-exp(-abs(lambda0)*theta0_2)-(lambda0*theta0_3));
matlabFunction(l, 'File', 'l_function');
%%

lambda0 = fsolve(@l_function,0);

%omega0 = (lambda0*v0 + v0)/R;
omega0 = - v0/(R*(lambda0-1));

tb0 = 2000;

theta0 = [theta0_1; theta0_2; theta0_3];

x0 = [v0; omega0; theta0];

u0 = tb0;

d0 = [wind0; slope0];

w0 = [d0; nu_w; nu_v];

y0 = [x0(2); x0(1)];

e0 = 0;

%% Linearized plant
%state
Amat = ABS_Amatrix(Cd,Ji,R,m,omega0,slope0,theta0_1,theta0_2,theta0_3,v0,wind0);
B1mat = ABS_B1matrix(Ji);
B2mat = ABS_B2matrix(Cd,Ji,R,m,omega0,slope0,theta0_1,theta0_2,theta0_3,v0,wind0);

%output
Cmat = ABS_Cmatrix(nu_w);
D1mat = ABS_D1matrix;
D2mat = ABS_D2matrix(omega0);

%error
Cemat = ABS_Cematrix(R,lambda_star,nu_w);
D1emat = ABS_D1ematrix;
D2emat = ABS_D2ematrix(R,lambda_star,omega0);

%% Jordan Canonical Form
[V,Vn,J] = JCF(Amat);

%% spider plot
n = length(Amat);
for kk = 1:n
    VV = Vn(kk,:);
    figure
    h = spider_plot(VV,...
        'AxesLimits', [zeros(1,n); 100*ones(1,n)],...
        'AxesLabels', {'$z_1$', '$z_2$','$z_3$','$z_4$','$z_5$'},...
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

rank(Rmat)
rank(Omat)

%% INITIAL CONDITIONS FOR THE NONLINEAR PLANT

v_init = v0; % [m/s] vehicle speed
omega_init = v_init/R; % [rad/s] rear wheel speed

x_init = [v_init; omega_init; theta0_1; theta0_2; theta0_3];

%% RUN THE SIMULATOR
PLANT = 0; % 0 = linear, 1 = nonlinear
TimeSpan = 200;
DT = 1e-6;
%% Simulink model
%out = sim('SimulinkModel',TimeSpan);
%save CurrentWorkspace

%% PLOT RESULTS

%run plotResults