function LHDrawTraj3D(rs, lambdaGs, Ls, tspan)
% Description:
% Create a 3D Plot of the propagated orbit in the ECEF reference frame.

global thetaG0

M = length(tspan);
rMatrixECI = zeros(M, 3);

D_sid = 86164;              % s
omegaE = 2*pi/D_sid;        % rad/s

for j = 1 : M

    r_vect = [rs(j), 0, 0];
    lambdaG = lambdaGs(j);
    L = Ls(j);

    t = tspan(j);
    thetaG = deg2rad(thetaG0) + omegaE * t;

    r_vect = r_vect * R(2, rad2deg(-L)) * R(3, rad2deg(lambdaG)) * R(3, rad2deg(thetaG));

    rMatrixECI(j, :) = r_vect;
    
end

X = rMatrixECI(:, 1);
Y = rMatrixECI(:, 2);
Z = rMatrixECI(:, 3);

[x,y,z]=sphere;

I = imread('earth.jpg');

earth = surface(6378.1363*x, 6378.1363*y, 6378.1363*z, flipud(I), 'FaceColor', 'texturemap', 'EdgeColor', 'none', 'CDataMapping', 'direct', 'EdgeAlpha', 1, 'FaceAlpha', 1);
earth_axis = [0 0 1];
rotate(earth, earth_axis, thetaG0)

hold on
plot3(X,Y,Z,'Color','#ff7403', 'Linestyle', '-', 'linewidth', 1.8)
plot3(X(1), Y(1), Z(1), 'Color', '#ff2e2e', 'LineStyle','none', 'marker', '.', 'markersize', 15)
plot3(X(end), Y(end), Z(end), 'Color', '#20facb', 'LineStyle','none', 'marker', '.', 'markersize', 15)
plot3(0,0,0,'g*')
hold off
grid on
axis equal
xlabel('x')
ylabel('y')
zlabel('z')
view([200, 20])


end