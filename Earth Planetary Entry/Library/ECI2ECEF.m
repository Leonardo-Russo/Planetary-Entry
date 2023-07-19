function [rMatrixECEF, thetaGs] = ECI2ECEF(rMatrixECI, tspan)
% Description: this function tranforms a set of ECI coordinates into ECEF
% from the known local parameters and the associated timespan vector.

global thetaG0

N = size(rMatrixECI, 1);
M = length(tspan);

if N ~= M
    error('Error! Input dimensions are NOT coherent!')
end

rMatrixECEF = zeros(M, 3);
thetaGs = zeros(M, 1);

D_sid = 86164;  % s
omega_E = rad2deg(2*pi/D_sid);   % deg/s

for i = 1 : M

    t = tspan(i);
    thetaG = thetaG0 + omega_E * t;
    
    rMatrixECEF(i, :) = rMatrixECI(i, :) * R(3, thetaG)';
    thetaGs(i) = thetaG;

end


end