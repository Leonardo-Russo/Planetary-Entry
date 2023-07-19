function dx = ECIDynamicalModel(t, x)
% Description: this function contains the Dynamical Model used for state
% integration between t0 and tEI.

N = length(x);
dx = zeros(N, 1);

mu = 398600.4415;   % km^3/s^2
r = norm(x(1:3));

dx(1:3) = x(4:6);
dx(4:6) = -mu/r^3*x(1:3);


end