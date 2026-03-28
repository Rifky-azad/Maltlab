function [u, du, strain, stress] = analytical_soln(x, p)
% Analytical solution for the spinning rod problem
%
% Governing ODE:  EA u''(x) + rho*A*omega^2*x = 0
% BC1: u(0) = 0  (fixed at center of rotation)
% BC2: EA*u'(L) = m*omega^2*L  (tip mass traction)
%
% Solution:
%   u(x) = (omega^2/E) * [ m*L*x/A + rho*L^2*x/2 - rho*x^3/6 ]
%
% Inputs:
%   x - vector of positions [m]
%   p - parameter struct from problem_params()
%
% Outputs:
%   u      - displacement [m]
%   du     - du/dx (strain = du/dx for 1D)
%   strain - engineering strain
%   stress - axial stress [Pa]

    w2 = p.omega^2;
    E  = p.E;
    rho = p.rho;
    L  = p.L;
    m  = p.m;
    A  = p.A;

    % Displacement
    u = (w2/E) * (m*L*x/A + rho*L^2*x/2 - rho*x.^3/6);

    % Strain = du/dx
    du = (w2/E) * (m*L/A + rho*L^2/2 - rho*x.^2/2);
    strain = du;

    % Stress = E * strain
    stress = E * strain;
end
