function plot_vector(O, v, color)
% Description: this function serves as a shortcut to plot a quiver3 vector
% given the origin, vector and color choices.

quiver3(O(1), O(2), O(3), v(1), v(2), v(3), 'Color', color)


end