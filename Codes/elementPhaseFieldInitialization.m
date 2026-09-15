function Kphi = ...
    elementPhaseFieldInitialization(mesh, e, material, quad)
% ------------------------------------------------------------
% elementPhaseFieldInitialization - Compute local AT2 initialization matrix
%
% PHASE-FIELD PHYSICS RELEVANCE:
% This function solves the element-level zero-load equilibrium problem:
%       phi - ell^2 * Laplacian(phi) = 0
%
% This PDE arises directly from minimizing the AT2 crack-surface functional 
% in the complete absence of elastic strain energy. By executing this over 
% the domain, the sharp discrete initial crack is mathematically smeared 
% into a diffuse damage band of width proportional to ell (l0).
%
% The weak form dictates that the local element matrix is:
%       Kphi = M_e + ell^2 * K_e
%
% where:
%       M_e = integral(N * N' dOmega)    -> Element Mass Matrix
%       K_e = integral(B * B' dOmega)    -> Element Diffusion Matrix
%
% Because this is a purely geometric/diffusive initialization, 
% no displacement field, history field, or mechanical degradation 
% function is required.
% ------------------------------------------------------------

    % --------------------------------------------------------
    % Element connectivity and coordinates
    % --------------------------------------------------------
    conn   = mesh.elem.conn(e,:);
    coords = mesh.nodes.coord(conn,:); % Physical coordinates of the element's nodes
    nen    = mesh.elem.nen;
    W      = quad.W;
    Q      = quad.Q;
    ngp    = numel(W);
    ell    = material.l0;
    
    % --------------------------------------------------------
    % Local matrix allocation
    % --------------------------------------------------------
    Kphi = zeros(nen,nen);
    
    % --------------------------------------------------------
    % Gauss-point numerical integration loop
    % --------------------------------------------------------
    for q = 1:ngp
        % Evaluate Shape functions (N) and their derivatives with respect 
        % to the reference/parent coordinates (dNdxi) at the current Gauss point.
        [N,dNdxi] = ...
            lagrange_basis(mesh.meta.elemType,Q(q,:));
            
        % ----------------------------------------------------
        % Isoparametric Mapping (Jacobian)
        % ----------------------------------------------------
        % J maps the derivatives from the reference domain (xi, eta) 
        % to the physical domain (x, y). 
        J = coords' * dNdxi;
        detJ = det(J);
        
        % Check Jacobian for element inversion/severe distortion.
        % A negative determinant implies the element geometry has folded 
        % in on itself, which will catastrophically corrupt the integration.
        if detJ <= 0
            error( ...
                'Non-positive Jacobian in element %d, GP %d.', ...
                e,q);
        end
        
        % Physical spatial derivatives (B-matrix for the scalar field)
        % dNdx = dNdxi * J^{-1}
        dNdx = dNdxi / J;
        
        % ----------------------------------------------------
        % Initialization operator assembly
        %
        % M_e = (N * N') * W_q * detJ
        % K_e = (dNdx * dNdx') * W_q * detJ
        %
        % Kphi = M_e + ell^2 * K_e
        % ----------------------------------------------------
        Kphi = Kphi + ...
            ( ...
                N*N' + ...
                ell^2*(dNdx*dNdx') ...
            ) ...
            * W(q) * detJ;
    end
end