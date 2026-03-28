function [L2_disp, L2_stress, energy_err] = compute_L2_error(result, p)
% Compute L2 error norms for displacement and stress, and energy norm error
%
% L2_disp   = sqrt( integral( (u_fem - u_exact)^2 dx ) )
% L2_stress = sqrt( integral( (sigma_fem - sigma_exact)^2 dx ) )
% energy_err = sqrt( integral( EA*(du_fem/dx - du_exact/dx)^2 dx ) )

    nElem = result.nElem;
    order = result.order;
    conn  = result.conn;
    x_nodes = result.x_nodes;
    U = result.U;

    nGauss = order + 2;  % use higher-order quadrature for error
    [gp, gw] = gauss_quad_local(nGauss);

    L2_disp_sq    = 0;
    L2_stress_sq  = 0;
    energy_err_sq = 0;

    for e = 1:nElem
        enodes = conn(e,:);
        xe = x_nodes(enodes);
        ue = U(enodes);

        for g = 1:nGauss
            xi = gp(g);
            w  = gw(g);

            [N, dNdxi] = shape_functions_local(xi, order);
            J = dNdxi' * xe;
            dNdx = dNdxi / J;

            x_gp = N' * xe;
            u_fem = N' * ue;
            du_fem = dNdx' * ue;

            [u_ex, du_ex, ~, sig_ex] = analytical_soln(x_gp, p);
            sig_fem = p.E * du_fem;

            L2_disp_sq    = L2_disp_sq   + (u_fem - u_ex)^2 * J * w;
            L2_stress_sq  = L2_stress_sq  + (sig_fem - sig_ex)^2 * J * w;
            energy_err_sq = energy_err_sq + p.EA * (du_fem - du_ex)^2 * J * w;
        end
    end

    L2_disp    = sqrt(L2_disp_sq);
    L2_stress  = sqrt(L2_stress_sq);
    energy_err = sqrt(energy_err_sq);
end


function [gp, gw] = gauss_quad_local(n)
    switch n
        case 1, gp=0; gw=2;
        case 2, gp=[-1;1]/sqrt(3); gw=[1;1];
        case 3, gp=[-sqrt(3/5);0;sqrt(3/5)]; gw=[5/9;8/9;5/9];
        case 4
            gp=[-sqrt(3/7+2/7*sqrt(6/5));-sqrt(3/7-2/7*sqrt(6/5));
                 sqrt(3/7-2/7*sqrt(6/5)); sqrt(3/7+2/7*sqrt(6/5))];
            gw=[(18-sqrt(30))/36;(18+sqrt(30))/36;
                (18+sqrt(30))/36;(18-sqrt(30))/36];
        case 5
            gp=[-1/3*sqrt(5+2*sqrt(10/7));-1/3*sqrt(5-2*sqrt(10/7));0;
                 1/3*sqrt(5-2*sqrt(10/7)); 1/3*sqrt(5+2*sqrt(10/7))];
            gw=[(322-13*sqrt(70))/900;(322+13*sqrt(70))/900;128/225;
                (322+13*sqrt(70))/900;(322-13*sqrt(70))/900];
        otherwise, error('Gauss quadrature not implemented for n=%d', n);
    end
end


function [N, dNdxi] = shape_functions_local(xi, order)
    switch order
        case 1
            N = [(1-xi)/2; (1+xi)/2];
            dNdxi = [-1/2; 1/2];
        case 2
            N = [xi*(xi-1)/2; 1-xi^2; xi*(xi+1)/2];
            dNdxi = [xi-1/2; -2*xi; xi+1/2];
        case 3
            nodes = [-1, -1/3, 1/3, 1];
            N = zeros(4,1); dNdxi = zeros(4,1);
            for i = 1:4
                den = 1; N_val = 1; dN_val = 0;
                for j = 1:4
                    if j ~= i, den = den*(nodes(i)-nodes(j)); end
                end
                for j = 1:4
                    if j ~= i
                        prod_others = 1;
                        for k = 1:4
                            if k~=i && k~=j
                                prod_others = prod_others*(xi-nodes(k));
                            end
                        end
                        dN_val = dN_val + prod_others;
                        N_val = N_val*(xi-nodes(j));
                    end
                end
                N(i) = N_val/den; dNdxi(i) = dN_val/den;
            end
        otherwise, error('Order %d not implemented', order);
    end
end
