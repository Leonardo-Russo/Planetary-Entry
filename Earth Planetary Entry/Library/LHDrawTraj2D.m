function LHDrawTraj3D(rs, lambdaGs, Ls, tspan)
% Description:
% Create a 3D Plot of the propagated orbit in the ECEF reference frame.

global thetaG0

M = length(tspan);
E = zeros(M, 1);

rE = 6378.136;      % km
h = rs - rE;        % km

for j = 1 : M-1

    if j < M
        r_avg = mean(rs(j:j+1));
        dlambdaG = lambdaGs(j+1) - lambdaGs(j);
        dE = r_avg * dlambdaG;
        E(j+1) = E(j) + dE;

    end

end


plot(E, h, 'LineStyle','-', 'LineWidth', 1.3, 'Color', '#ff7403')
xlabel('East Direction $\ [km]$', 'Interpreter','latex', 'FontSize', 12)
ylabel('Altitude $\ [km]$', 'Interpreter','latex', 'FontSize', 12)


end