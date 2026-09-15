function psiPlusGP = ...
    elementPositiveEnergy( ...
    mesh,e,u,material,model,quad)
% ============================================================
% elementPositiveEnergy - Extract Active Tensile Energy
%
% PHASE-FIELD HPC RELEVANCE:
% After the displacement Newton-Raphson solver converges, the framework 
% must evaluate the instantaneous active tensile strain energy (psi_plus) 
% across the entire domain to update the thermodynamic history field (H).
%
% While the elementDisplacementContribution routine also evaluates the 
% split schemes, it carries the massive computational overhead of 
% assembling the 4th-order consistent tangent stiffness matrix (K_e) and 
% the internal residual force vector (R_e). 
%
% To maximize computational efficiency, this dedicated subroutine acts as 
% a "read-only" function. It computes the local kinematics, extracts the 
% scalar active strain energy from the split scheme, and entirely bypasses 
% all matrix assembly operations.
% ============================================================

dim = mesh.meta.dim;

% ------------------------------------------------------------
% Element data extraction
% ------------------------------------------------------------
conn   = mesh.elem.conn(e,:);
coords = mesh.nodes.coord(conn,:);
dofU   = mesh.elem.dofU(e,:);
u_e    = u(dofU);
nen    = mesh.elem.nen;

W = quad.W;
Q = quad.Q;
ngp = numel(W);

% Pre-allocate the vector to store the scalar energy at each Gauss point
psiPlusGP = zeros(ngp,1);

%% ============================================================
% Gauss integration loop
% ============================================================
for q = 1:ngp
    
    %% --------------------------------------------------------
    % Shape functions and Reference Derivatives
    % ---------------------------------------------------------
    [~,dNdxi] = ...
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
    
    % Physical spatial derivatives
    dNdx = dNdxi / J;
    
    %% --------------------------------------------------------
    % Strain-Displacement (B) Matrix
    % ---------------------------------------------------------
    B = displacementBMatrix(dNdx,dim,nen);
    
    %% --------------------------------------------------------
    % Current Mechanical Strain
    % ---------------------------------------------------------
    strain = B*u_e;
    
    %% --------------------------------------------------------
    % Tension-Compression Constitutive Split
    % ---------------------------------------------------------
    % This calls the core constitutive router to decompose the intact 
    % strain energy based on the user's selected split scheme 
    % ('none', 'amor', or 'spectral').
    split = computeSplitInfo(strain,material,model);
    
    %% --------------------------------------------------------
    % Positive Energy Extraction
    % ---------------------------------------------------------
    % We strictly extract the active/tensile component of the strain 
    % energy density. The compressive/passive component is discarded 
    % here, as it does not contribute to crack propagation.
    psiPlusGP(q) = split.psiPlus;
end
end