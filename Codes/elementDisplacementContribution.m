function [Ke,Re] = ...
    elementDisplacementContribution( ...
    mesh,e,u,phi,material,model,quad)
% ============================================================
% elementDisplacementContribution - Mechanical Tangent and Residual
%
% PHASE-FIELD PHYSICS RELEVANCE:
% This function evaluates the weak form of the mechanical equilibrium 
% equation at the element level. Because the solver relies on a full 
% Newton-Raphson scheme for the displacement field, this routine computes 
% both the internal force vector (Re) and the consistent tangent stiffness 
% matrix (Ke) simultaneously.
%
% HPC RELEVANCE:
% This function is called millions of times during a full crack propagation 
% analysis (at every element, for every Newton iteration, for every load step). 
% Therefore, it is stripped of all non-essential calculations: 
% no energy calculations, no stress storage, and no history tracking. 
% It strictly computes what is necessary for the Newton system assembly.
% ============================================================

dim = mesh.meta.dim;

%% ------------------------------------------------------------
% Element data extraction
% ------------------------------------------------------------
conn   = mesh.elem.conn(e,:);
coords = mesh.nodes.coord(conn,:);
dofU   = mesh.elem.dofU(e,:);

% Extract local nodal variables from the global state vectors
u_e   = u(dofU);
phi_e = phi(conn);
nen   = mesh.elem.nen;
ndofu = dim*nen;

W = quad.W;
Q = quad.Q;
ngp = numel(W);

%% ------------------------------------------------------------
% Initialize local matrices
% ------------------------------------------------------------
Ke = zeros(ndofu,ndofu);
Re = zeros(ndofu,1);

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
    
    % Physical spatial derivatives of shape functions
    dNdx = dNdxi / J;
    
    %% --------------------------------------------------------
    % Strain-Displacement (B) Matrix
    % ---------------------------------------------------------
    % Constructs the standard small-strain symmetric gradient operator.
    B = displacementBMatrix(dNdx,dim,nen);
    
    %% --------------------------------------------------------
    % Current Mechanical Strain
    % ---------------------------------------------------------
    % Evaluates the infinitesimal strain tensor (in Voigt form) at the GP.
    strain = B*u_e;
    
    %% --------------------------------------------------------
    % Current Phase Field (Damage)
    % ---------------------------------------------------------
    % Interpolates the nodal damage values to the Gauss point.
    phi_gp = N'*phi_e;
    
    %% --------------------------------------------------------
    % Degradation Function
    % ---------------------------------------------------------
    % g(phi) reduces the stiffness of the material as phi -> 1. 
    % A residual stiffness k is inherently included inside this function 
    % to prevent det(Ctan) = 0 when fully fractured.
    [g,~,~] = ...
        degradationFunction(phi_gp,model,material);
        
    %% --------------------------------------------------------
    % Constitutive Split (Tension-Compression Asymmetry)
    % ---------------------------------------------------------
    % Computes the active/passive stress tensors (sigma) and 
    % tangent moduli (C) based on the user's chosen split scheme.
    split = computeSplitInfo(strain,material,model);
    
    %% --------------------------------------------------------
    % Stress and Consistent Tangent Assembly
    % ---------------------------------------------------------
    % To physically prevent crack faces from interpenetrating under compression, 
    % the degradation function g(phi) is multiplied ONLY by the tensile/active 
    % portions of the stress and stiffness tensors. The compressive portions 
    % remain pristine and undegraded.
    switch lower(model.split)
        case 'none'
            % Isotropic: Degrade everything equally.
            sigma = g*split.sigma;
            Ctan  = g*split.C;
            
        case {'amor','spectral'}
            % Asymmetric: Degrade only the positive (tensile/volumetric) parts.
            sigma = g*split.sigmaPlus + split.sigmaMinus;
            Ctan  = g*split.Cplus + split.Cminus;
            
        otherwise
            error('Unknown split: %s', model.split);
    end
    
    %% --------------------------------------------------------
    % Convert stress tensor to FEM Voigt form
    % ---------------------------------------------------------
    % Maps the 2x2 or 3x3 Cauchy stress tensor into a column vector 
    % compatible with the B-matrix for the dot product.
    sigmaVoigt = stressTensorToVoigt(sigma,dim);
    
    %% --------------------------------------------------------
    % Integration weight (dV)
    % ---------------------------------------------------------
    weight = W(q)*detJ;
    
    %% --------------------------------------------------------
    % Internal residual vector (Out-of-balance forces)
    % Re = Integral( B^T * sigma * dV )
    % ---------------------------------------------------------
    Re = Re + B'*sigmaVoigt*weight;
    
    %% --------------------------------------------------------
    % Consistent Newton Tangent Stiffness Matrix
    % Ke = Integral( B^T * Ctan * B * dV )
    % ---------------------------------------------------------
    Ke = Ke + B'*Ctan*B*weight;
end
end