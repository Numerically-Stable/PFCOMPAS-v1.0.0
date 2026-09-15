function [Eel,Efrac] = ...
    elementEnergy( ...
    mesh,e,u,phi, ...
    material,model,quad)
% ============================================================
% elementEnergy - Local Integration of System Energies
%
% PHASE-FIELD PHYSICS RELEVANCE:
% This subroutine numerically integrates the two fundamental components 
% of the total potential energy functional over a single finite element:
%
%   1. Eel (Elastic Energy): The mechanical strain energy stored in the 
%      bulk material. It accounts for tension-compression asymmetry, 
%      meaning only the active (tensile) portion is degraded by damage.
%   2. Efrac (Fracture Energy): The energy dissipated to create the 
%      diffuse crack surface, governed by the AT2 regularization profile.
%
% HPC RELEVANCE:
% This is a "read-only" post-processing function. It strictly avoids 
% constructing the B-matrix into a global stiffness array, focusing 
% exclusively on evaluating the scalar energy densities at the Gauss points 
% and accumulating them over the element volume.
% ============================================================

conn   = mesh.elem.conn(e,:);
coords = mesh.nodes.coord(conn,:);
nen    = mesh.elem.nen;
dim    = model.dim;

u_e   = u(mesh.elem.dofU(e,:));
phi_e = phi(conn);

% Initialize local element energy scalars
Eel   = 0;
Efrac = 0;

%% ============================================================
% Gauss Integration Loop
% ============================================================
for q = 1:numel(quad.W)
    
    %% --------------------------------------------------------
    % Shape functions and Isoparametric Mapping
    % ---------------------------------------------------------
    [N,dNdxi] = ...
        lagrange_basis( ...
        mesh.meta.elemType, ...
        quad.Q(q,:));
        
    J = coords' * dNdxi;
    detJ = det(J);
    dNdx = dNdxi / J;
    
    %% --------------------------------------------------------
    % B matrix and Mechanical Strain
    % ---------------------------------------------------------
    B = displacementBMatrix(dNdx,dim,nen);
    strain = B*u_e;
    
    %% --------------------------------------------------------
    % Phase Field and Spatial Gradient
    % ---------------------------------------------------------
    phiGP = N' * phi_e;
    
    % gradPhi is a column vector [dphi/dx; dphi/dy; (dphi/dz)]
    gradPhi = dNdx' * phi_e;
    
    %% --------------------------------------------------------
    % Tension-Compression Split (Energy Decomposition)
    % ---------------------------------------------------------
    split = computeSplitInfo(strain, material, model);
    
    %% --------------------------------------------------------
    % Degradation Function
    % ---------------------------------------------------------
    [g,~,~] = degradationFunction(phiGP,model,material);
    
    %% --------------------------------------------------------
    % Degraded Elastic Strain Energy Density (psi_e)
    % ---------------------------------------------------------
    % Irrespective of whether Amor or Spectral splits are used, the 
    % degradation g(phi) is applied ONLY to the active tensile energy.
    % The compressive/passive energy (psiMinus) is completely preserved.
    psi = g*split.psiPlus + split.psiMinus;
    
    %% --------------------------------------------------------
    % Regularized Fracture Surface Density (psi_frac)
    % ---------------------------------------------------------
    % For the AT2 model (normalization constant c0 = 2):
    % psi_frac = (Gc / (2*l0)) * phi^2 + (Gc*l0 / 2) * |grad(phi)|^2
    %
    % Note: gradPhi'*gradPhi is the exact dot product yielding 
    % the squared L2 norm of the phase-field spatial gradient.
    Gc = material.Gc;
    l0 = material.l0;
    
    psiFrac = ...
        Gc/(2*l0)*phiGP^2 ...
        + ...
        Gc*l0/2*(gradPhi'*gradPhi);
        
    %% --------------------------------------------------------
    % Volumetric Integration (dV)
    % ---------------------------------------------------------
    weight = quad.W(q)*detJ;
    
    % Accumulate element energies
    Eel   = Eel + psi*weight;
    Efrac = Efrac + psiFrac*weight;
end
end