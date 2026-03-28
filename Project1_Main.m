%% ========================================================================
%  MAE 404/503 — Project 1: FEA of a Spinning Aluminum Rod
%  ========================================================================
%  An aluminum rod (d=5cm, E=70GPa, rho=2700 kg/m^3) swings a 1 kg mass
%  around a circle of radius 1 m at 2000 RPM.
%
%  Body force:  b(x) = rho*A*omega^2*x
%  BC: u(0) = 0,  EA*u'(L) = m*omega^2*L
%  ========================================================================

clear; close all; clc;
p = problem_params();

fprintf('=== Project 1: Spinning Rod FEA ===\n');
fprintf('E = %.0f GPa,  rho = %.0f kg/m^3\n', p.E/1e9, p.rho);
fprintf('d = %.0f mm,   A = %.4e m^2\n', p.d*1000, p.A);
fprintf('L = %.1f m,    m = %.1f kg\n', p.L, p.m);
fprintf('omega = %.2f rad/s (%.0f RPM)\n', p.omega, p.RPM);
fprintf('Tip force F = %.2f N\n\n', p.F_tip);

% Fine grid for analytical solution
x_fine = linspace(0, p.L, 500)';
[u_exact, ~, ~, sig_exact] = analytical_soln(x_fine, p);

%% ====== C1: Problem Statement and Assumptions ======
fprintf('--- C1: Problem Statement ---\n');
fprintf('Geometry: 1D rod, x in [0, L], L = 1 m\n');
fprintf('Material: Aluminum, E = 70 GPa, rho = 2700 kg/m^3\n');
fprintf('Cross-section: circular, d = 5 cm, A = %.4e m^2\n', p.A);
fprintf('Loading: centripetal body force b(x) = rho*A*omega^2*x\n');
fprintf('         Tip traction: F = m*omega^2*L = %.2f N\n', p.F_tip);
fprintf('BCs: u(0) = 0 (fixed),  EA*u''(L) = F (natural)\n');
fprintf('Assumptions:\n');
fprintf('  - Small deformations (linear elasticity)\n');
fprintf('  - Uniform cross-section\n');
fprintf('  - 1D axial deformation only\n');
fprintf('  - Steady-state rotation (constant omega)\n\n');

%% ====== C2: Analytical Displacement ======
fprintf('--- C2: Analytical Displacement ---\n');
fprintf('u(x) = (omega^2/E)*[m*L*x/A + rho*L^2*x/2 - rho*x^3/6]\n');
fprintf('u(L) = %.6e m\n\n', u_exact(end));

%% ====== C3: Analytical Strain and Stress ======
fprintf('--- C3: Analytical Strain and Stress ---\n');
fprintf('strain(x) = du/dx = (omega^2/E)*[m*L/A + rho*L^2/2 - rho*x^2/2]\n');
fprintf('stress(x) = E*strain(x) = omega^2*[m*L/A + rho*L^2/2 - rho*x^2/2]\n');
fprintf('stress(0) = %.4f MPa\n', sig_exact(1)/1e6);
fprintf('stress(L) = %.4f MPa\n\n', sig_exact(end)/1e6);

%% ====== C4: Linear FEM (2 & 4 elements) — Displacement ======
fprintf('--- C4: Linear FEM Displacement ---\n');
res_lin2 = fem_1d(2, 1, p);
res_lin4 = fem_1d(4, 1, p);
fprintf('Linear 2-elem: u(L) = %.6e m\n', res_lin2.U(end));
fprintf('Linear 4-elem: u(L) = %.6e m\n\n', res_lin4.U(end));

figure('Name','C4: Linear FEM Displacement');
plot(x_fine, u_exact, 'k-', 'LineWidth', 1.5); hold on;
plot(res_lin2.x_nodes, res_lin2.U, 'bo-', 'LineWidth', 1.2, 'MarkerSize', 8);
plot(res_lin4.x_nodes, res_lin4.U, 'rs-', 'LineWidth', 1.2, 'MarkerSize', 6);
xlabel('x [m]'); ylabel('u(x) [m]');
title('C4: Linear FEM Displacement (2 & 4 Elements)');
legend('Analytical', 'Linear 2-elem', 'Linear 4-elem', 'Location', 'northwest');
grid on;

%% ====== C5: Linear FEM Stress ======
fprintf('--- C5: Linear FEM Stress ---\n');
[x_s2, sig_s2, ~] = fem_stress(res_lin2, p);
[x_s4, sig_s4, ~] = fem_stress(res_lin4, p);

figure('Name','C5: Linear FEM Stress');
plot(x_fine, sig_exact/1e6, 'k-', 'LineWidth', 1.5); hold on;
plot(x_s2, sig_s2/1e6, 'b-', 'LineWidth', 1.2);
plot(x_s4, sig_s4/1e6, 'r-', 'LineWidth', 1.2);
xlabel('x [m]'); ylabel('\sigma(x) [MPa]');
title('C5: Linear FEM Stress (Piecewise Constant per Element)');
legend('Analytical', 'Linear 2-elem', 'Linear 4-elem', 'Location', 'northeast');
grid on;

%% ====== C6: Quadratic FEM (2 & 4 elements) — Displacement ======
fprintf('--- C6: Quadratic FEM Displacement ---\n');
res_quad2 = fem_1d(2, 2, p);
res_quad4 = fem_1d(4, 2, p);
fprintf('Quadratic 2-elem: u(L) = %.6e m\n', res_quad2.U(end));
fprintf('Quadratic 4-elem: u(L) = %.6e m\n\n', res_quad4.U(end));

figure('Name','C6: Quadratic FEM Displacement');
plot(x_fine, u_exact, 'k-', 'LineWidth', 1.5); hold on;
% For quadratic, interpolate through nodes for smooth curve
for e = 1:res_quad2.nElem
    enodes = res_quad2.conn(e,:);
    xe = res_quad2.x_nodes(enodes);
    ue = res_quad2.U(enodes);
    xi_plot = linspace(-1,1,50);
    x_plot = zeros(50,1); u_plot = zeros(50,1);
    for i = 1:50
        N = [xi_plot(i)*(xi_plot(i)-1)/2; 1-xi_plot(i)^2; xi_plot(i)*(xi_plot(i)+1)/2];
        x_plot(i) = N'*xe; u_plot(i) = N'*ue;
    end
    if e == 1
        plot(x_plot, u_plot, 'b-', 'LineWidth', 1.2);
    else
        plot(x_plot, u_plot, 'b-', 'LineWidth', 1.2, 'HandleVisibility','off');
    end
end
plot(res_quad2.x_nodes, res_quad2.U, 'bo', 'MarkerSize', 8);

for e = 1:res_quad4.nElem
    enodes = res_quad4.conn(e,:);
    xe = res_quad4.x_nodes(enodes);
    ue = res_quad4.U(enodes);
    xi_plot = linspace(-1,1,50);
    x_plot = zeros(50,1); u_plot = zeros(50,1);
    for i = 1:50
        N = [xi_plot(i)*(xi_plot(i)-1)/2; 1-xi_plot(i)^2; xi_plot(i)*(xi_plot(i)+1)/2];
        x_plot(i) = N'*xe; u_plot(i) = N'*ue;
    end
    if e == 1
        plot(x_plot, u_plot, 'r-', 'LineWidth', 1.2);
    else
        plot(x_plot, u_plot, 'r-', 'LineWidth', 1.2, 'HandleVisibility','off');
    end
end
plot(res_quad4.x_nodes, res_quad4.U, 'rs', 'MarkerSize', 6);

xlabel('x [m]'); ylabel('u(x) [m]');
title('C6: Quadratic FEM Displacement (2 & 4 Elements)');
legend('Analytical', 'Quad 2-elem', 'Quad 2 nodes', 'Quad 4-elem', 'Quad 4 nodes', 'Location', 'northwest');
grid on;

%% ====== C7: Quadratic FEM Stress ======
fprintf('--- C7: Quadratic FEM Stress ---\n');
[x_sq2, sig_sq2, ~] = fem_stress(res_quad2, p);
[x_sq4, sig_sq4, ~] = fem_stress(res_quad4, p);

figure('Name','C7: Quadratic FEM Stress');
plot(x_fine, sig_exact/1e6, 'k-', 'LineWidth', 1.5); hold on;
plot(x_sq2, sig_sq2/1e6, 'b-', 'LineWidth', 1.2);
plot(x_sq4, sig_sq4/1e6, 'r-', 'LineWidth', 1.2);
xlabel('x [m]'); ylabel('\sigma(x) [MPa]');
title('C7: Quadratic FEM Stress (Correct Trend)');
legend('Analytical', 'Quad 2-elem', 'Quad 4-elem', 'Location', 'northeast');
grid on;

%% ====== C8: Overlay — Analytical vs FEM ======
fprintf('--- C8: Overlay Analytical vs FEM ---\n');
res_lin_fine = fem_1d(20, 1, p);
res_quad_fine = fem_1d(20, 2, p);

figure('Name','C8: Analytical vs FEM Overlay');
plot(x_fine, u_exact, 'k-', 'LineWidth', 2); hold on;
plot(res_lin_fine.x_nodes, res_lin_fine.U, 'b--o', 'LineWidth', 1, 'MarkerSize', 4);
plot(res_quad_fine.x_nodes, res_quad_fine.U, 'r--s', 'LineWidth', 1, 'MarkerSize', 4);
xlabel('x [m]'); ylabel('u(x) [m]');
title('C8: Analytical vs FEM Displacement (20 elements)');
legend('Analytical', 'Linear FEM (20 elem)', 'Quadratic FEM (20 elem)', 'Location', 'northwest');
grid on;

%% ====== V1: Mesh Refinement Study ======
fprintf('\n--- V1: Mesh Refinement Study ---\n');
nElems_v1 = [2, 4, 8, 16, 32, 64, 128];
fprintf('%-10s %-15s %-15s %-15s %-15s\n', 'nElem', 'max |u| Lin', 'max |u| Quad', 'max|sig| Lin', 'max|sig| Quad');

for i = 1:length(nElems_v1)
    ne = nElems_v1(i);
    rl = fem_1d(ne, 1, p);
    rq = fem_1d(ne, 2, p);
    [~, sl, ~] = fem_stress(rl, p);
    [~, sq, ~] = fem_stress(rq, p);
    fprintf('%-10d %-15.6e %-15.6e %-15.4f %-15.4f\n', ne, max(rl.U), max(rq.U), max(sl)/1e6, max(sq)/1e6);
end
fprintf('\n');

%% ====== V2: L2 Error vs Mesh Size (Linear) ======
fprintf('--- V2: L2 Error vs Mesh Size (Linear) ---\n');
nElems_conv = [2, 4, 8, 16, 32, 64, 128, 256];
h_lin = zeros(size(nElems_conv));
L2_lin = zeros(size(nElems_conv));
energy_lin = zeros(size(nElems_conv));

for i = 1:length(nElems_conv)
    res = fem_1d(nElems_conv(i), 1, p);
    [L2_lin(i), ~, energy_lin(i)] = compute_L2_error(res, p);
    h_lin(i) = res.h;
end

% Estimate convergence rate
coeffs_lin = polyfit(log(h_lin), log(L2_lin), 1);
rate_lin = coeffs_lin(1);
fprintf('Linear elements: L2 displacement convergence rate ~ %.2f\n', rate_lin);

figure('Name','V2: L2 Error — Linear');
loglog(h_lin, L2_lin, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 8); hold on;
loglog(h_lin, h_lin.^2 * L2_lin(1)/h_lin(1)^2, 'k--', 'LineWidth', 1);
xlabel('Element size h [m]'); ylabel('L_2 displacement error');
title('V2: L2 Error vs Mesh Size (Linear Elements)');
legend(sprintf('Linear FEM (slope=%.2f)', rate_lin), 'O(h^2) reference', 'Location', 'southeast');
grid on;

%% ====== V3: L2 Error vs Mesh Size (Quadratic) ======
fprintf('--- V3: L2 Error vs Mesh Size (Quadratic) ---\n');
h_quad = zeros(size(nElems_conv));
L2_quad = zeros(size(nElems_conv));
energy_quad = zeros(size(nElems_conv));

for i = 1:length(nElems_conv)
    res = fem_1d(nElems_conv(i), 2, p);
    [L2_quad(i), ~, energy_quad(i)] = compute_L2_error(res, p);
    h_quad(i) = res.h;
end

coeffs_quad = polyfit(log(h_quad), log(L2_quad), 1);
rate_quad = coeffs_quad(1);
fprintf('Quadratic elements: L2 displacement convergence rate ~ %.2f\n\n', rate_quad);

figure('Name','V3: L2 Error — Quadratic');
loglog(h_quad, L2_quad, 'rs-', 'LineWidth', 1.5, 'MarkerSize', 8); hold on;
loglog(h_quad, h_quad.^3 * L2_quad(1)/h_quad(1)^3, 'k--', 'LineWidth', 1);
xlabel('Element size h [m]'); ylabel('L_2 displacement error');
title('V3: L2 Error vs Mesh Size (Quadratic Elements)');
legend(sprintf('Quadratic FEM (slope=%.2f)', rate_quad), 'O(h^3) reference', 'Location', 'southeast');
grid on;

%% ====== Combined Convergence Plot ======
figure('Name','Convergence Comparison');
loglog(h_lin, L2_lin, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 8); hold on;
loglog(h_quad, L2_quad, 'rs-', 'LineWidth', 1.5, 'MarkerSize', 8);
loglog(h_lin, energy_lin, 'b^--', 'LineWidth', 1.2, 'MarkerSize', 6);
loglog(h_quad, energy_quad, 'rv--', 'LineWidth', 1.2, 'MarkerSize', 6);
xlabel('Element size h [m]'); ylabel('Error norm');
title('Convergence: Linear vs Quadratic Elements');
legend(sprintf('L2 Linear (%.2f)', rate_lin), sprintf('L2 Quad (%.2f)', rate_quad), ...
       'Energy Linear', 'Energy Quad', 'Location', 'southeast');
grid on;

fprintf('=== All tasks complete. ===\n');
