function spin_globe()
% Description: this function serves no purpose other than creating a
% spinning globe with a given time step provided in seconds.


[x,y,z]=sphere;

I = imread('earth.jpg');

earth = surface(6378.1363*x, 6378.1363*y, 6378.1363*z, flipud(I), 'FaceColor', 'texturemap', 'EdgeColor', 'none', 'CDataMapping', 'direct', 'EdgeAlpha', 1, 'FaceAlpha', 1);
earth_axis = [0 0 1];


hold on
plot3(0,0,0,'g*')
hold off
grid on
axis equal
xlabel('x')
ylabel('y')
zlabel('z')
view([120, 20])

theta = 0;

while true
    
    theta = theta + 0.001;
    rotate(earth, earth_axis, theta)
    pause(0.01)

end

end