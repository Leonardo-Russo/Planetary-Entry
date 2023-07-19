function ECEFDrawTraj3D(rMatrixECI, tspan)
% Description:
% Create a 3D Plot of the propagated orbit in the ECI reference frame.

rMatrixECEF = ECI2ECEF(rMatrixECI, tspan);

X = rMatrixECEF(:, 1);
Y = rMatrixECEF(:, 2);
Z = rMatrixECEF(:, 3);

[x,y,z]=sphere;

I = imread('earth.jpg');

surface(6378.1363*x, 6378.1363*y, 6378.1363*z, flipud(I), 'FaceColor', 'texturemap', 'EdgeColor', 'none', 'CDataMapping', 'direct', 'EdgeAlpha', 0.75, 'FaceAlpha', 0.75)

hold on
plot3(X,Y,Z,'Color','#ff7403', 'Marker', '.','MarkerSize', 5, 'Linestyle', 'none')
plot3(0,0,0,'g*')
hold off
grid on
axis equal
xlabel('x')
ylabel('y')
zlabel('z')
view([120, 20])


end