function [Ls, lambdaGs] = ECEF2LH(rMatrixECEF)
% Description: this function converts from ECEF to LH frame.

M = size(rMatrixECEF, 1);

Ls = zeros(M, 1);
lambdaGs = zeros(M, 1);

for i = 1 : M
    
    r_vect = rMatrixECEF(i, :);
    x = r_vect(1);
    y = r_vect(2);
    z = r_vect(3);
    r = norm(r_vect);

    Ls(i) = rad2deg(asin(z/r));
    lambdaGs(i) = rad2deg(atan2(y, x));

end

end