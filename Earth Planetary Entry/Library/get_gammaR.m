function [gammaR, lambdaG] = get_gammaR(r, v, t)
% Description: this function retrieves the flight path angle gammaR from
% the provided input.

[r_vect_ECEF, thetaG] = ECI2ECEF(r', t);
[L, lambdaG] = ECEF2LH(r_vect_ECEF);

v_vectECEF = v' * R(3, thetaG)';
v_vectLH = v_vectECEF * R(3, lambdaG)' * R(2, -L)';

vr = v_vectLH(1);
vt = v_vectLH(2);

gammaR = rad2deg(atan2(vr, vt));


end