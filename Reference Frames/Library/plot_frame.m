function plot_frame(F, O, color)
% Description: this function plots the reference frames from their base in
% F and origin point in O.

plot3(O(1), O(2), O(3), 'Color','c', 'marker', '*')
hold on

quiver3(O(1), O(2), O(3), F(1, 1), F(2, 1), F(3, 1), "Color", color, 'LineWidth', 1.5)
quiver3(O(1), O(2), O(3), F(1, 2), F(2, 2), F(3, 2), "Color", color, 'LineWidth', 1.5)
quiver3(O(1), O(2), O(3), F(1, 3), F(2, 3), F(3, 3), "Color", color, 'LineWidth', 1.5)

text(O(1) + F(1, 1), O(2) + F(2, 1), O(3) + F(3, 1), '1', 'Color', color)
text(O(1) + F(1, 2), O(2) + F(2, 2), O(3) + F(3, 2), '2', 'Color', color)
text(O(1) + F(1, 3), O(2) + F(2, 3), O(3) + F(3, 3), '3', 'Color', color)

end