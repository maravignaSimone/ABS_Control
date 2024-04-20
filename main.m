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

slope0 = 0;
wind0 = 0;
nu_w = 0;
nu_v = 0;
v0 = 28; % [m/s] initial velocity

% friction coefficients
%dry asphalt
theta0_1=1.28;
theta0_2=23.99;
theta0_3=0.52;
lambda_star = 0.17;
%wet asphalt
%theta0_1 = 0.86;
%theta0_2 = 33.82;
%theta0_3 = 0.35;
%lambda_star=0.131;

% function for finding lambda0 of the equilibrium point
l = (Cd*(v0-wind0)^2) + m*g*sign(lambda0)*theta0_1*(1-exp(-abs(lambda0)*theta0_2)-lambda0*theta0_3);
matlabFunction(l, 'File', 'l_function');
%%

lam = fsolve(@l_function,0);

omega0 = (lam*v0 + v0)/R;

%%

Amat = ABS_Amatrix(Cd,Ji,R,m,omega0,slope0,theta0_1,theta0_2,theta0_3,v0,wind0);
B1mat = ABS_B1matrix(Ji);
B2mat = ABS_B2matrix(Cd,Ji,R,m,omega0,slope0,theta0_1,theta0_2,theta0_3,v0,wind0);
Cmat = ABS_Cmatrix(nu_w);
D1mat = ABS_D1matrix;
D2mat = ABS_D2matrix(omega0);
Cemat = ABS_Cematrix(R,lambda_star,nu_w);
D1emat = ABS_D1ematrix;
D2emat = ABS_D2ematrix(R,lambda_star,omega0);

%%
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

%%
Rmat = ctrb(Amat,B1mat);
Omat = obsv(Amat,Cmat);

rank(Rmat)
rank(Omat)

%%
PLANT = 0;
tilde_x_init = [0; 0; 0; 0; 0];
u0 = 0;
w0 = [0 0 0 0]';

%% RUN THE SIMULATOR
TimeSpan = 100;
DT = 1e-8;
out = sim('SimulinkModel',TimeSpan);
save CurrentWorkspace

%% PLOT RESULTS

run plotResults