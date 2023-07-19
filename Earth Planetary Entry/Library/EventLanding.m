function [value, isterminal, direction] = EventLanding(t, x)
% Description: this function is the event function which detects when the
% SC lands on the Earth Surface.

rE = 6378.136;          % km

r = x(1);         % km

value = r - rE;         % value that needs to be null
isterminal = 1;         % stop the integration?
direction = 0;          % reach value from both directions?


end