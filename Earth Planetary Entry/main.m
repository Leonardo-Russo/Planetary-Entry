%% Homework 6 - Leonardo Russo 2015563

close all
clear all
clc

addpath('Library/')
addpath('Output/')

%% Options

dv_case = '1';

savechoice = '0';
additional_plots = '1';

% Define the Options for ode113()
Tol0 = 1e-9;
Tol1 = 1e-11;
MaxStep = 1e-1;

optionsODE_ECI = odeset('Events', @EventEI, 'RelTol', Tol0, 'AbsTol',Tol1);

optionsODE_Reentry = odeset('Events', @EventLanding, 'RelTol', Tol0, 'AbsTol',Tol1, 'MaxStep', MaxStep);

%% Define Known Quantities

global thetaG0
thetaG0 = 60;

m = 2400;       % kg
S = 3.8;        % m^2
Cd = 1.341;     % avg value
L_tilde = 0;    % Lift
Rc = 2.235;     % m

h0 = 400;       % km
rE = 6378.136;      % km
mu = 398600.4415;   % km^3/s^2
D_sid = 86164;              % s
D_sid = 86400;
omegaE = 2*pi/D_sid;        % rad/s
omegaE_vect = [0, 0, omegaE]';

t0 = 0;             % s
tEI_max = 100000;   % s

h_EI = 100;     % km


% Define Initial State of the S/C
if dv_case == '1'
    dv = 0.2;   % km/s
elseif dv_case == '2'
    dv = 0.4;   % km/s
end

lambdaG0 = 30;      % deg

r0_vectECI = R(3, thetaG0 + lambdaG0)' * [rE+h0, 0, 0]';
r0_ECI = norm(r0_vectECI);

v0_vectECI = R(3, thetaG0 + lambdaG0)' * [0, sqrt(mu/r0_ECI) - dv, 0]';


%% Part A - ECI

% Define Initial Conditions
x0 = [r0_vectECI; v0_vectECI];

% State Propagation up to EI Conditions
[tspan, x, te, xe, ~] = ode113(@ECIDynamicalModel, [t0, tEI_max], x0, optionsODE_ECI);

% Store Propagation Results
M = length(tspan);
rMatrix = x(:, 1:3);
vMatrix = x(:, 4:6);

rI_ini = xe(1:3)';
vI_ini = xe(4:6)';
tEI = tspan(end);

gammaR = get_gammaR(rI_ini, vI_ini, tEI);    % deg


%%% Results of Part A %%%

figure('Name', 'Initial Descent of S/C towards EI')

ECIDrawTraj3D(rMatrix)
if savechoice == '1'
    saveas(gcf, strcat('Output\ECIprop-', dv_case,'.jpg'))
end

var_namesA = ["vI,ini", "gammaI,ini"]';
unitsA = ["km/s", "deg"]';
summaryA = [norm(vI_ini), gammaR]';

fprintf('\n\n\t\t<strong>Summary of Part A</strong>\n\n')
disp(table(var_namesA, unitsA, summaryA, 'VariableNames',["Quantity", "Units", "Value"]))


%% Part B

vR_ini = vI_ini - cross(omegaE_vect, rI_ini);

[gammaR, lambdaG] = get_gammaR(rI_ini, vR_ini, tEI);        % deg


%%% Results of Part B %%%

var_namesB = ["vR,ini", "gammaR,ini"]';
unitsB = ["km/s", "deg"]';
summaryB = [norm(vR_ini), gammaR]';

fprintf('\n\n\t\t<strong>Summary of Part B</strong>\n\n')
disp(table(var_namesB, unitsB, summaryB, 'VariableNames',["Quantity", "Units", "Value"]))


%% Part C

r_EI = norm(rI_ini);
vR_EI = norm(vR_ini);
lambdaG_EI = deg2rad(lambdaG);
gammaR_AE = deg2rad(gammaR);

% Define Initial Conditions
x0 = [r_EI, lambdaG_EI, gammaR_AE, vR_EI];
t0 = tEI;                   % s
tLand_max = te + 10000;     % s

% State Propagation up to Landing
[tspan, x, te, xe, ie] = ode113(@B1DynamicalModel, [t0, tLand_max], x0, optionsODE_Reentry);


% Store Propagation Results
M = length(tspan);
r = x(:, 1);
lambdaG = x(:, 2);
L = zeros(M, 1);
gammaR = x(:, 3);
vR = x(:, 4);
hi = r - rE;
tLand = te;


%%% Results of Part C

figure('Name', 'Reentry of the S/C using B1 Dynamical Model')

LHDrawTraj3D(r, lambdaG, L, tspan)
if savechoice == '1'
    saveas(gcf, strcat('Output\B1reentry-', dv_case,'.jpg'))
end


figure('Name', 'Altitude Evolution in Time')

plot(tspan, hi, 'LineStyle','-', 'LineWidth', 1.3, 'Color', '#ff7403')
xlabel('$t \ [s]$', 'Interpreter','latex', 'FontSize', 12)
ylabel('$h \ [km]$', 'Interpreter','latex', 'FontSize', 12)
if savechoice == '1'
    saveas(gcf, strcat('Output\h_t-', dv_case,'.jpg'))
end


figure('Name', 'Relative Flight Path Angle Evolution in Time')

plot(tspan, rad2deg(gammaR), 'LineStyle','-', 'LineWidth', 1.3, 'Color', '#ff7403')
xlabel('$t \ [s]$', 'Interpreter','latex', 'FontSize', 12)
ylabel('$\gamma_R \ [deg]$', 'Interpreter','latex', 'FontSize', 12)
if savechoice == '1'
    saveas(gcf, strcat('Output\gammaR_t-', dv_case,'.jpg'))
end


figure('Name', 'Relative Velocity Evolution in Time')

plot(tspan, vR, 'LineStyle','-', 'LineWidth', 1.3, 'Color', '#ff7403')
xlabel('$t \ [s]$', 'Interpreter','latex', 'FontSize', 12)
ylabel('$v_R \ [km/s]$', 'Interpreter','latex', 'FontSize', 12)
if savechoice == '1'
    saveas(gcf, strcat('Output\vR_t-', dv_case,'.jpg'))
end


figure('Name', 'Reentry Trajectory')

LHDrawTraj2D(r, lambdaG, L, tspan)
if savechoice == '1'
    saveas(gcf, strcat('Output\reentry_trajectory-', dv_case,'.jpg'))
end


%% Part D


R_A = [0 1 0;...
       0 0 -1;...
       -1 0 0];

M = length(tspan);
q0 = zeros(M, 1);
q = zeros(M, 3);

yaw = zeros(M, 1);
pitch = zeros(M, 1);
roll = zeros(M, 1);

for i = 1 : M

    thetaGi = thetaG0 + rad2deg(omegaE) * tspan(i);
    lambdaGi = rad2deg(lambdaG(i));
    gammaRi = rad2deg(gammaR(i));

    RNC = R_A * R(3, -gammaRi) * R(3, lambdaGi) * R(3, thetaGi);

    [q0(i), q(i, :)] = C2q(RNC);

end

eulers = quat2eul([q0, q]);
psi = rad2deg(eulers(:, 1));
theta = rad2deg(eulers(:, 2));
phi = rad2deg(eulers(:, 3));


figure('Name', 'Quaternions Evolution in Time')

hold on
plot(tspan, q0, 'LineStyle','-', 'LineWidth', 1.3)
plot(tspan, q(:, 1), 'LineStyle','-', 'LineWidth', 1.3)
plot(tspan, q(:, 2), 'LineStyle','-', 'LineWidth', 1.3)
plot(tspan, q(:, 3), 'LineStyle','-', 'LineWidth', 1.3)
xlabel('$t \ [s]$', 'Interpreter','latex', 'FontSize', 12)
ylabel('$q_i$', 'Interpreter','latex', 'FontSize', 12)
legend('$q_0$', '$q_1$', '$q_2$', '$q_3$', 'interpreter', 'latex', 'fontsize', 12, 'location', 'best')
if savechoice == '1'
    saveas(gcf, strcat('Output\quaternions-', dv_case,'.jpg'))
end


if additional_plots == '1'
    figure('Name', 'Euler Angles Evolution in Time')
    
    hold on
    plot(tspan, psi, 'LineStyle','-', 'LineWidth', 1.3)
    plot(tspan, theta, 'LineStyle','-', 'LineWidth', 1.3)
    plot(tspan, phi, 'LineStyle','-', 'LineWidth', 1.3)
    xlabel('$t \ [s]$', 'Interpreter','latex', 'FontSize', 12)
    ylabel('Euler Angles $\ [deg]$', 'Interpreter','latex', 'FontSize', 12)
    legend('$\psi$ - Yaw', '$\theta$ - Pitch', '$\varphi$ - Roll', 'interpreter', 'latex', 'fontsize', 10, 'location', 'best')
    if savechoice == '1'
        saveas(gcf, strcat('Output\angles-', dv_case,'.jpg'))
    end
end


%% Part E

g0 = 9.8065 * 1e-3;     % km/s^2
K = 1.762e-4;           % kg^1/2 m^-1

TOF = tLand - tEI;

% Initialize Local Quantities
aD = zeros(M, 1);
Pd = zeros(M, 1);
qs = zeros(M, 1);

% Compute Local Quantities
for i = 1 : M

    ri = r(i);
    lambdaGi = lambdaG(i);
    gammaRi = gammaR(i);
    vRi = vR(i);

    rho0 = 1.225;               % kg/m^3
    Beta = 0.13718;            % km^-1
    hi = ri - rE;
    rho = rho0 * exp(-Beta*hi);
    D = 0.5*Cd*S*rho*vRi^2 * 1e3;

    ridot = vRi*sin(gammaRi);
    lambdaGidot = vRi*cos(gammaRi) / ri;
    gammaRidot = (-mu/(ri^2*vRi) + vRi/ri)*cos(gammaRi) + omegaE^2*ri/vRi*cos(gammaRi) + 2*omegaE;
    vRidot = -mu/ri^2*sin(gammaRi) - D/m + omegaE^2*ri*sin(gammaRi);

    aD(i) = -0.5 / sqrt(vRi^2 + omegaE^2*ri^2 + 2*omegaE*vRi*ri*cos(gammaRi)) * (...
            2*vRi*vRidot + 2*ri*omegaE^2*ridot + 2*omegaE*(ri*cos(gammaRi)*vRidot +...
            vRi*cos(gammaRi)*ridot - ri*vRi*sin(gammaRi)*gammaRidot) );  % km/s^2


    Pd(i) = 0.5 * rho * vRi^2 * 1e1;   % bar

    qs(i) = K/sqrt(Rc) * sqrt(rho) * vRi^3 * 1e3;     % MW/m^2

end

[aDmax, iaDmax] = max(aD);
[Pdmax, iPdmax] = max(Pd);
[qsmax, iqsmax] = max(qs);
h = r - rE;


%%% Results of Part E %%%

var_namesE = ["Time of Flight", "aD,max", "Pd,max", "qs,max"]';
unitsE = ["s", "g0", "bar", "MW/m^2"]';
summaryE = [TOF, aDmax/g0, Pdmax, qsmax]';
timesE = ["---", tspan(iaDmax)-tEI, tspan(iPdmax)-tEI, tspan(iqsmax)-tEI]';

fprintf('\n\n\t\t<strong>Summary of Part E</strong>\n\n')
disp(table(var_namesE, unitsE, summaryE, timesE, 'VariableNames',["Quantity", "Units", "Value", "Time"]))



%% Part F

% Introduce the Constants for the Allen-Eggers Solution
gammaR_AE = gammaR(1);
BC = Cd*S/m;                % m^2/kg
C = BC*rho0/(2*Beta*sin(gammaR_AE)) * 1e3;

% Compute the Analytical Allen-Eggers Solution
vR_AE = zeros(M, 1);
qs_AE = zeros(M, 1);
aD_AE = zeros(M, 1);
Pd_AE = zeros(M, 1);

for i = 1 : M

    hi = h(i);
    rho = rho0 * exp(-Beta*hi);
    vR_AEi = vR_EI * exp(C * exp(-Beta*hi));

    vR_AE(i) = vR_AEi;      % km/s
    qs_AE(i) = K/sqrt(Rc) * sqrt(rho) * vR_AEi^3 * 1e3;     % MW/m^2
    aD_AE(i) = BC*rho0/2 * vR_EI^2 * exp(-Beta*hi) * exp(2*C*exp(-Beta*hi)) * 1e3;  % km/s^2
    Pd_AE(i) = 0.5 * rho * vR_AEi^2 * 1e1;   % bar

end


%%% Results of Part F

figure('Name', 'Relative Flight Path Angle Evolution with Altitude')

hold on
plot(h, rad2deg(gammaR_AE*ones(M, 1)), 'LineStyle','-', 'LineWidth', 1.3, 'Color', '#39c7fa')
plot(h, rad2deg(gammaR), 'LineStyle','-', 'LineWidth', 1.3, 'Color', '#ff7403')
xlabel('$h \ [km]$', 'Interpreter','latex', 'FontSize', 12)
ylabel('$\gamma_R \ [deg]$', 'Interpreter','latex', 'FontSize', 12)
legend('Allen-Eggers Solution', 'Analytical Solution', 'fontsize', 10, 'location', 'best')
if savechoice == '1'
    saveas(gcf, strcat('Output\gammaR_h-', dv_case,'.jpg'))
end


figure('Name', 'Relative Velocity Evolution with Altitude')

hold on
plot(h, vR_AE, 'LineStyle','-', 'LineWidth', 1.3, 'Color', '#39c7fa')
plot(h, vR, 'LineStyle','-', 'LineWidth', 1.3, 'Color', '#ff7403')
xlabel('$h \ [km]$', 'Interpreter','latex', 'FontSize', 12)
ylabel('$v_R \ [km/s]$', 'Interpreter','latex', 'FontSize', 12)
legend('Allen-Eggers Solution', 'Analytical Solution', 'fontsize', 10, 'location', 'best')
if savechoice == '1'
    saveas(gcf, strcat('Output\vR_h-', dv_case,'.jpg'))
end


figure('Name', 'Stagnation Point Intensive Heat Evolution with Altitude')

hold on
plot(h, qs_AE, 'LineStyle','-', 'LineWidth', 1.3, 'Color', '#39c7fa')
plot(h, qs, 'LineStyle','-', 'LineWidth', 1.3, 'Color', '#ff7403')
xlabel('$h \ [km]$', 'Interpreter','latex', 'FontSize', 12)
ylabel('$q_s \ [MW/m^2]$', 'Interpreter','latex', 'FontSize', 12)
legend('Allen-Eggers Solution', 'Analytical Solution', 'fontsize', 10, 'location', 'best')
if savechoice == '1'
    saveas(gcf, strcat('Output\qs_h-', dv_case,'.jpg'))
end


figure('Name', 'Inertial Deceleration Evolution with Altitude')

hold on
plot(h, aD_AE*1e3, 'LineStyle','-', 'LineWidth', 1.3, 'Color', '#39c7fa')
plot(h, aD*1e3, 'LineStyle','-', 'LineWidth', 1.3, 'Color', '#ff7403')
xlabel('$h \ [km]$', 'Interpreter','latex', 'FontSize', 12)
ylabel('$a_D \ [m/s^2]$', 'Interpreter','latex', 'FontSize', 12)
legend('Allen-Eggers Solution', 'Analytical Solution', 'fontsize', 10, 'location', 'best')
if savechoice == '1'
    saveas(gcf, strcat('Output\aD_h-', dv_case,'.jpg'))
end


figure('Name', 'Dynamical Pressure Evolution with Altitude')

hold on
plot(h, Pd_AE, 'LineStyle','-', 'LineWidth', 1.3, 'Color', '#39c7fa')
plot(h, Pd, 'LineStyle','-', 'LineWidth', 1.3, 'Color', '#ff7403')
xlabel('$h \ [km]$', 'Interpreter','latex', 'FontSize', 12)
ylabel('$P_d \ [bar]$', 'Interpreter','latex', 'FontSize', 12)
legend('Allen-Eggers Solution', 'Analytical Solution', 'fontsize', 10, 'location', 'best')
if savechoice == '1'
    saveas(gcf, strcat('Output\Pd_h-', dv_case,'.jpg'))
end

