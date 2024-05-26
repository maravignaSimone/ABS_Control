%Plot of the friction coefficient in dry asphalt changing Theta1 parameter of the 10%
clc
clear
close all

lambda = linspace(-1, 1, 200);
i = 1;
% friction coeffiients
% 1) Dry asphalt
% 2) wet asphalt
% 3) Snow
% 4) Ice
% 5) Dry Cobblestone
% 6) wet cobblestone
theta1 = [1.28 0.86 0.19 1.37 0.4];
theta2 = [23.99 33.82 94.13 6.46 33.71];
theta3 = [0.52 0.35 0.05 0.67 0.12];

mean_theta1 = mean(theta1);
mean_theta2 = mean(theta2);
mean_theta3 = mean(theta3);
std_theta1 = std(theta1);
std_theta2 = std(theta2);
std_theta3 = std(theta3);

% mean_theta1 = theta1(i);
% mean_theta2 = theta2(i);
% mean_theta3 = theta3(i);

hold on
mu_l = sign(lambda) .* mean_theta1 .* (1 - exp(-abs(lambda) .* mean_theta2)) - lambda .* mean_theta3;
mu_l_thetaplus10 = sign(lambda) .* (mean_theta1+std_theta1) .* (1 - exp(-abs(lambda) .* mean_theta2)) - lambda .* mean_theta3;
mu_l_thetaminus10 = sign(lambda) .* (mean_theta1-std_theta1) .* (1 - exp(-abs(lambda) .* mean_theta2)) - lambda .* mean_theta3;
% mu_l_thetaplus10 = sign(lambda) .* (mean_theta1) .* (1 - exp(-abs(lambda) .* (mean_theta2+std_theta2))) - lambda .* mean_theta3;
% mu_l_thetaminus10 = sign(lambda) .* (mean_theta1) .* (1 - exp(-abs(lambda) .* (mean_theta2-std_theta2))) - lambda .* mean_theta3;
% mu_l_thetaplus10 = sign(lambda) .* (mean_theta1) .* (1 - exp(-abs(lambda) .* mean_theta2)) - lambda .* (mean_theta3-std_theta3);
% mu_l_thetaminus10 = sign(lambda) .* (mean_theta1) .* (1 - exp(-abs(lambda) .* mean_theta2)) - lambda .* (mean_theta3-std_theta3);

plot(lambda, mu_l, 'LineWidth', 1)
plot(lambda, mu_l_thetaplus10, 'LineWidth', 1)
plot(lambda, mu_l_thetaminus10, 'LineWidth', 1)

xlabel('\lambda')
ylabel('\mu_l(\lambda, \Theta)')
title('Plot of \mu(\lambda, \Theta)')
show