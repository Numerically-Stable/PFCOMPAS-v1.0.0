function [Kphi,Rphi] = ...
    elementPhaseFieldContribution( ...
    mesh,e,phi,H,model,material,quad)
% ============================================================
% elementPhaseFieldContribution - Phase-Field Tangent and Residual
%
% PHASE-FIELD PHYSICS RELEVANCE:
% This subroutine evaluates the local weak form of the phase-field 
% evolution equation (analogous to an Allen-Cahn or Ginzburg-Landau PDE).
% 
% The residual Rphi represents the energetic competition between:
%   1. The thermodynamic driving force trying to increase damage 
%      (driven by the strain history field H and g'(phi)).
%   2. The fracture resistance trying to prevent damage 
%      (driven by the material toughness Gc and regularized by l0).
%
% HPC RELEVANCE:
% This function strictly computes the arrays necessary for the Newton 
% solver. It bypasses any global fracture energy (dissipation) calculations, 
% deferring those non-essential post-processing steps until after global 
% equilibrium is completely satisfied.
% ============================================================

conn   = mesh.elem.conn(e,:);
coords = mesh.nodes.coord(conn,:);
phi_e  = phi(conn);

% Extract the history field strictly at the Gauss points for this element.
% Avoiding nodal interpolation of H preserves the sharpness of the crack tip.
H_e    = H(e,:);

nen = mesh.elem.nen;
W   = quad.W;
Q   = quad.Q;
ngp = numel(W);

Gc  = material.Gc;
ell = material.l0;

% ------------------------------------------------------------
% AT2 Normalization Constant
% In the Francfort-Marigo / Bourdin variational approach, the fracture 
% surface density functional is normalized by a constant c0 to ensure that 
% the spatial integral of the diffuse crack exactly converges to Gc as l0 -> 0.
% For the AT2 model (quadratic geometric function), c0 = 2.
% For the AT1 model (linear geometric function), c0 = 8/3.
% ------------------------------------------------------------
c0 = 2;

%% ============================================================
% Initialize local matrices
% ============================================================
Kphi = zeros(nen,nen);
Rphi = zeros(nen,1);

%% ============================================================
% Gauss integration loop
% ============================================================
for q = 1:ngp
    %% --------------------------------------------------------
    % Shape functions and Reference Derivatives
    % ---------------------------------------------------------
    [N,dNdxi] = ...
        lagrange_basis( ...
        mesh.meta.elemType,Q(q,:));
        
    %% --------------------------------------------------------
    % Geometry and Isoparametric Jacobian
    % ---------------------------------------------------------
    J = coords' * dNdxi;
    detJ = det(J);
    if detJ <= 0
        error('Non-positive Jacobian: element %d, GP %d.', e,q);
    end
    dNdx = dNdxi / J;
    
    %% --------------------------------------------------------
    % Phase field and its spatial gradient
    % ---------------------------------------------------------
    phi_gp = N'*phi_e;
    gradPhi = dNdx'*phi_e;
    
    %% --------------------------------------------------------
    % History field at current integration point
    % ---------------------------------------------------------
    H_gp = H_e(q);
    
    %% --------------------------------------------------------
    % Degradation derivatives
    % ---------------------------------------------------------
    % Retrieves g'(phi) and g''(phi). 
    % Note: g(phi) itself is not needed here; it only scales the mechanics.
    [~,dg,ddg] = ...
        degradationFunction(phi_gp,model,material);
        
    %% --------------------------------------------------------
    % AT2 Crack Geometric Function w(phi)
    %
    % For AT2: w(phi) = phi^2
    % dw  = d(w)/d(phi) = 2*phi
    % ddw = d2(w)/d(phi)^2 = 2
    % ---------------------------------------------------------
    dw  = 2*phi_gp;
    ddw = 2;
    
    %% --------------------------------------------------------
    % Integration weight (dV)
    % ---------------------------------------------------------
    weight = W(q)*detJ;
    
    %% --------------------------------------------------------
    % Internal residual (Out-of-balance energy forces)
    %
    % Rphi = Integral [ N * { H * g'(phi) + (Gc/(c0*l0)) * w'(phi) } 
    %                   + Grad(N) * { (2*Gc*l0/c0) * Grad(phi) } ] dV
    % ---------------------------------------------------------
    Rphi = Rphi ...
        + ...
        ( ...
        N*( ...
            H_gp*dg ...
            + Gc/(c0*ell)*dw ...
        ) ...
        + ...
        dNdx*gradPhi*(2*Gc*ell/c0) ...
        )*weight;
        
    %% --------------------------------------------------------
    % Consistent Newton Tangent (Jacobian K_phi)
    %
    % Kphi = Integral [ N * N^T * { H * g''(phi) + (Gc/(c0*l0)) * w''(phi) }
    %                   + Grad(N) * Grad(N)^T * { (2*Gc*l0/c0) } ] dV
    % ---------------------------------------------------------
    Kphi = Kphi ...
        + ...
        ( ...
        N*N'*( ...
            H_gp*ddg ...
            + Gc/(c0*ell)*ddw ...
        ) ...
        + ...
        dNdx*dNdx'*(2*Gc*ell/c0) ...
        )*weight;
end
end