function p = problem_params()
% Returns a struct with all physical parameters for Project 1
%
% Aluminum rod swinging a 1 kg mass at 2000 RPM in a circle of radius 1 m

    p.E     = 70e9;          % Young's modulus [Pa]
    p.rho   = 2700;          % Density [kg/m^3]
    p.d     = 0.05;          % Diameter [m]
    p.A     = pi*(p.d/2)^2;  % Cross-sectional area [m^2]
    p.L     = 1.0;           % Rod length (radius of circle) [m]
    p.m     = 1.0;           % Tip mass [kg]
    p.RPM   = 2000;          % Rotational speed [RPM]
    p.omega = p.RPM * 2*pi / 60;  % Angular velocity [rad/s]

    % Derived quantities
    p.EA = p.E * p.A;
    p.F_tip = p.m * p.omega^2 * p.L;  % Traction at x=L from point mass [N]
end
