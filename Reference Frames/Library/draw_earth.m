function draw_earth()
% Description: this function draws the earth and sets it as background.

plot3(0, 0, 0, 'Color','c', 'MarkerSize',0.1)

axis equal

[x,y,z]=sphere;

I = imread('earth.jpg');

surface(6378.1363*x, 6378.1363*y, 6378.1363*z, flipud(I), 'FaceColor', 'texturemap', 'EdgeColor', 'none', 'CDataMapping', 'direct', 'EdgeAlpha', 0.8, 'FaceAlpha', 0.8)

view([120, 25])


end