%% Initialize Workspace

close all
clear all
clc

addpath('Library/')

rE = 6371;      % km
axis_length = 2000;       % km

O = [0, 0, 0];

%% Define the Frames

ECI = 2*rE * eye(3);

theta_G = 20;   % Greenwhich Latitude
ECEF = ECI * R(3, theta_G)';

L = 20;         % Latitude
lambda_G = 10;      % Geographical Longitude

SC = [8000, 0, 0] * R(2, -L) * R(3, lambda_G) * R(3, theta_G);

LH = axis_length/(2*rE) * ECEF * R(3, lambda_G)' * R(2, -L)';

xi_R = 40;      % Heading Angle
AO = LH * R(1, xi_R)';

gamma_R = 30;   % Flight Path Angle
RV = AO * R(3, -gamma_R)';

sigma = 30;     % Bank Angle
WA = RV * R(2, sigma)';

alpha = 60;     % Angle of Attack
beta = 30;      % Sideslip Angle
ABA = WA * R(1, beta)' * R(3, -alpha)';

R_A = [0 1 0;...
       0 0 -1;...
       -1 0 0];
BA = ABA * R_A;


%% Plot the Frames

draw_earth()
hold on

plot_frame(ECI, O, 'b')
plot_frame(ECEF, O, 'g')
plot_frame(LH, SC, 'r')
plot_plane(LH(:, 2), LH(:, 3), SC)
plot_frame(RV, SC, 'k')
plot_frame(ABA, SC, 'c')



