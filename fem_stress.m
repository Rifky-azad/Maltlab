function [x_stress, stress_vals, strain_vals] = fem_stress(result, p)
% Compute stress and strain from FEM solution
%
% For each element, evaluates stress at Gauss points (or uniformly spaced points)
% Returns piecewise stress for plotting
%
% Inputs:
%   result - output from fem_1d()
%   p      - parameter struct
%
% Outputs:
%   x_stress   - x-coordinates for stress evaluation
%   stress_vals - stress values [Pa]
%   strain_vals - strain values

    nElem = result.nElem;
    order = result.order;
    conn  = result.conn;
    x_nodes = result.x_nodes;
    U = result.U;

    nPts = 20;  % evaluation points per element
    x_stress = zeros(nElem * nPts, 1);
    stress_vals = zeros(nElem * nPts, 1);
    strain_vals = zeros(nElem * nPts, 1);

    for e = 1:nElem
        enodes = conn(e,:);
        xe = x_nodes(enodes);
        ue = U(enodes);

        xi_eval = linspace(-1, 1, nPts)';
        idx = (e-1)*nPts + (1:nPts);

        for i = 1:nPts
            xi = xi_eval(i);
            [N, dNdxi] = shape_functions_local(xi, order);

            J = dNdxi' * xe;
            dNdx = dNdxi / J;

            x_stress(idx(i)) = N' * xe;
            strain_vals(idx(i)) = dNdx' * ue;
            stress_vals(idx(i)) = p.E * strain_vals(idx(i));
        end
    end
end


function [N, dNdxi] = shape_functions_local(xi, order)
% Same shape functions as in fem_1d (duplicated to keep files self-contained)
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
        otherwise
            error('Order %d not implemented', order);
    end
end
