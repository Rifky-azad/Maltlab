function result = fem_1d(nElem, order, p, varargin)
% 1D Finite Element solver for the spinning rod problem
%
% Inputs:
%   nElem  - number of elements
%   order  - polynomial order (1=linear, 2=quadratic, 3=cubic)
%   p      - parameter struct from problem_params()
%   Optional name-value pairs:
%     'sparse', true/false  - use sparse matrices (default: false)
%     'nGauss', n           - number of Gauss points (default: order+1)
%
% Outputs:
%   result.x_nodes  - global node coordinates
%   result.U        - nodal displacements
%   result.K        - global stiffness matrix
%   result.F        - global force vector
%   result.conn     - element connectivity
%   result.nDOF     - total degrees of freedom
%   result.h        - element size

    % Parse optional arguments
    opts = struct('sparse', false, 'nGauss', order+1);
    for i = 1:2:length(varargin)
        opts.(varargin{i}) = varargin{i+1};
    end

    nGauss = opts.nGauss;
    L = p.L;
    EA = p.EA;
    rho = p.rho;
    A = p.A;
    omega = p.omega;

    % Generate mesh
    h = L / nElem;
    nNodesPerElem = order + 1;
    nNodes = nElem * order + 1;

    % Global node coordinates (equally spaced)
    x_nodes = linspace(0, L, nNodes)';

    % Element connectivity (1-based indexing)
    conn = zeros(nElem, nNodesPerElem);
    for e = 1:nElem
        conn(e,:) = (e-1)*order + (1:nNodesPerElem);
    end

    % Gauss quadrature on [-1,1]
    [gp, gw] = gauss_quad(nGauss);

    % Assemble global stiffness and force
    nDOF = nNodes;
    if opts.sparse
        K = sparse(nDOF, nDOF);
    else
        K = zeros(nDOF, nDOF);
    end
    F = zeros(nDOF, 1);

    for e = 1:nElem
        enodes = conn(e,:);
        xe = x_nodes(enodes);  % element node coordinates

        Ke = zeros(nNodesPerElem);
        Fe = zeros(nNodesPerElem, 1);

        for g = 1:nGauss
            xi = gp(g);
            w  = gw(g);

            % Shape functions and derivatives in parent coords
            [N, dNdxi] = shape_functions(xi, order);

            % Jacobian: dx/dxi
            J = dNdxi' * xe;

            % Shape function derivatives in physical coords
            dNdx = dNdxi / J;

            % Physical coordinate at this Gauss point
            x_gp = N' * xe;

            % Body force at Gauss point: b(x) = rho*A*omega^2*x
            b_gp = rho * A * omega^2 * x_gp;

            % Element stiffness
            Ke = Ke + (EA * (dNdx * dNdx') * J * w);

            % Element force
            Fe = Fe + (N * b_gp * J * w);
        end

        % Assemble
        K(enodes, enodes) = K(enodes, enodes) + Ke;
        F(enodes) = F(enodes) + Fe;
    end

    % Apply tip load at x=L (last node)
    F(nDOF) = F(nDOF) + p.F_tip;

    % Apply essential BC: u(0) = 0 (node 1)
    % Penalty method or direct elimination — use elimination
    K_free = K(2:end, 2:end);
    F_free = F(2:end);

    % Solve
    U_free = K_free \ F_free;
    U = [0; U_free];

    % Pack results
    result.x_nodes = x_nodes;
    result.U       = U;
    result.K       = K;
    result.F       = F;
    result.conn    = conn;
    result.nDOF    = nDOF;
    result.h       = h;
    result.order   = order;
    result.nElem   = nElem;
end


function [gp, gw] = gauss_quad(n)
% Returns Gauss-Legendre quadrature points and weights on [-1,1]
    switch n
        case 1
            gp = 0;
            gw = 2;
        case 2
            gp = [-1; 1]/sqrt(3);
            gw = [1; 1];
        case 3
            gp = [-sqrt(3/5); 0; sqrt(3/5)];
            gw = [5/9; 8/9; 5/9];
        case 4
            gp = [-sqrt(3/7 + 2/7*sqrt(6/5));
                   -sqrt(3/7 - 2/7*sqrt(6/5));
                    sqrt(3/7 - 2/7*sqrt(6/5));
                    sqrt(3/7 + 2/7*sqrt(6/5))];
            gw = [(18 - sqrt(30))/36;
                  (18 + sqrt(30))/36;
                  (18 + sqrt(30))/36;
                  (18 - sqrt(30))/36];
        case 5
            gp = [-1/3*sqrt(5 + 2*sqrt(10/7));
                  -1/3*sqrt(5 - 2*sqrt(10/7));
                   0;
                   1/3*sqrt(5 - 2*sqrt(10/7));
                   1/3*sqrt(5 + 2*sqrt(10/7))];
            gw = [(322 - 13*sqrt(70))/900;
                  (322 + 13*sqrt(70))/900;
                  128/225;
                  (322 + 13*sqrt(70))/900;
                  (322 - 13*sqrt(70))/900];
        otherwise
            error('Gauss quadrature not implemented for n=%d', n);
    end
end


function [N, dNdxi] = shape_functions(xi, order)
% Lagrange shape functions and derivatives on [-1,1]
%
% order=1: linear    (2 nodes at xi=-1,1)
% order=2: quadratic (3 nodes at xi=-1,0,1)
% order=3: cubic     (4 nodes at xi=-1,-1/3,1/3,1)

    switch order
        case 1
            N = [(1-xi)/2; (1+xi)/2];
            dNdxi = [-1/2; 1/2];

        case 2
            N = [xi*(xi-1)/2;
                 1 - xi^2;
                 xi*(xi+1)/2];
            dNdxi = [xi - 1/2;
                     -2*xi;
                     xi + 1/2];

        case 3
            % Nodes at xi = -1, -1/3, 1/3, 1
            xi1 = -1; xi2 = -1/3; xi3 = 1/3; xi4 = 1;
            N = zeros(4,1);
            dNdxi = zeros(4,1);
            nodes = [xi1, xi2, xi3, xi4];
            for i = 1:4
                num = 1; dnum = 0; den = 1;
                for j = 1:4
                    if j ~= i
                        den = den * (nodes(i) - nodes(j));
                    end
                end
                % Compute N_i and dN_i/dxi using product rule
                N_val = 1;
                dN_val = 0;
                for j = 1:4
                    if j ~= i
                        % Product rule for derivative
                        prod_others = 1;
                        for k = 1:4
                            if k ~= i && k ~= j
                                prod_others = prod_others * (xi - nodes(k));
                            end
                        end
                        dN_val = dN_val + prod_others;
                        N_val = N_val * (xi - nodes(j));
                    end
                end
                N(i) = N_val / den;
                dNdxi(i) = dN_val / den;
            end

        otherwise
            error('Order %d not implemented', order);
    end
end
