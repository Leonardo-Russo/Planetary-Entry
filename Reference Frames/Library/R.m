function [RotMat] = R(axis, angle)
% Description: this function returns the rotation matrix around the input
% axis of input angle provided in degrees.

if axis == 1

    RotMat = [1, 0, 0;...
              0, cosd(angle), sind(angle);...
              0, -sind(angle), cosd(angle)];

elseif axis == 2

    RotMat = [cosd(angle), 0, -sind(angle);...
              0, 1, 0;...
              sind(angle), 0, cosd(angle)];

elseif axis == 3

    RotMat = [cosd(angle), sind(angle), 0;...
              -sind(angle), cosd(angle), 0;...
              0, 0, 1];

else

    error('Sorry, provided axis is not valid!\n')

end

end