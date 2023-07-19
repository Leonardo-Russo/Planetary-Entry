function dx = B1DynamicalModel(~, x)
% Description: this function hold the full dynamical model for the reentry
% problem.

% Import Data from Input
r = x(1);
lambdaG = x(2);
gammaR = x(3);
vR = x(4);

% Initialize State Update
N = length(x);
dx = zeros(N, 1);

% Define Known Quantities
m = 2400;                   % kg
Cd = 1.341;
S = 3.8;                    % m^2
L_tilde = 0;
sigma = 0;

rE = 6378.136;              % km
D_sid = 86164;              % s
omegaE = 2*pi/D_sid;        % rad/s
mu = 398600.4415;           % km^3/s^2

rho0 = 1.225;               % kg/m^3
B_rho = 0.13718;            % km^-1

% Compute Local Quantities
h = r - rE;
rho = rho0 * exp(-B_rho*h);
D = 0.5*Cd*S*rho*vR^2 * 1e3;


% Define State Update
dx(1) = vR*sin(gammaR);
dx(2) = vR*cos(gammaR) / r;
dx(3) = (-mu/(r^2*vR) + vR/r)*cos(gammaR) + (L_tilde*cos(sigma))/(m*vR) + omegaE^2*r*cos(gammaR)/vR + 2*omegaE;
dx(4) = -mu/r^2*sin(gammaR) - D/m + omegaE^2*r*sin(gammaR);


end