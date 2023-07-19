function [value, isterminal, direction] = EventEI(t, x)
% Description: this function is the event function which detects when the
% SC enters arrives at the EI.

rE = 6378.136;  % km
h_EI = 100;     % km

r_EI = rE + h_EI;       % km

r = norm(x(1:3));

value = r - r_EI;       % value that needs to be null
isterminal = 1;         % stop the integration?
direction = 0;          % reach value from both directions?


end