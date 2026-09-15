function [K,R] = ...
    assembleDisplacementSystem( ...
    mesh,u,phi,material,model,quad, ...
    parallel,assemblyPattern)
% ============================================================
% assembleDisplacementSystem - Global Mechanical Assembly
%
% PHASE-FIELD HPC RELEVANCE:
% This function drives the core computational effort of the mechanical 
% Newton-Raphson solver. Because the spectral or Amor split makes the 
% local stiffness matrix highly non-linear, this assembly is called 
% repeatedly (at every Newton iteration).
%
% The routine strictly enforces the parallel HPC paradigm:
% 1. PARFOR computes only the element-level physical integrations.
% 2. No threads write to global memory concurrently (avoiding race conditions).
% 3. The global sparse pattern (I, J, R indices) supplied from the 
%    preprocessing phase is used to perform a single-shot vectorized assembly.
% ============================================================

nelem = mesh.elem.n;
nen   = mesh.elem.nen;
ndofu = mesh.meta.dim*nen; % Degrees of freedom per element for displacement
ndof  = mesh.meta.dim*mesh.nodes.n; % Total global displacement DOFs

%% ============================================================
% Element numerical values (Pre-allocation)
% ============================================================
% These are the ONLY quantities accumulated by the parallel workers.
% V_K holds the flattened tangent stiffness matrices.
% V_R holds the local residual force vectors.
V_K = zeros(ndofu*ndofu,nelem);
V_R = zeros(ndofu,nelem);

%% ============================================================
% Element loop (Integration Phase)
% ============================================================
if parallel.useParfor
    parfor e = 1:nelem
        % Compute local tangent stiffness (Ke) and internal forces (Re)
        % based on the current displacement (u) and damage (phi) state.
        [Ke,Re] = ...
            elementDisplacementContribution( ...
            mesh,e,u,phi,material,model,quad);
            
        % Store values in independent column slices to guarantee thread safety
        V_K(:,e) = Ke(:);
        V_R(:,e) = Re;
    end
else
    for e = 1:nelem
        [Ke,Re] = ...
            elementDisplacementContribution( ...
            mesh,e,u,phi,material,model,quad);
        V_K(:,e) = Ke(:);
        V_R(:,e) = Re;
    end
end

%% ============================================================
% Global sparse assembly (Stiffness Matrix K)
% ============================================================
% Vectorized summation of overlapping node stiffness contributions.
K = sparse( ...
    assemblyPattern.u.I, ...
    assemblyPattern.u.J, ...
    V_K(:), ...
    ndof,ndof);

%% ============================================================
% Global residual (Internal Force Vector R)
% ============================================================
% While sparse() handles 2D matrices, accumarray() is the optimal HPC 
% equivalent for 1D vectors in MATLAB. It groups the flattened local 
% residual values V_R(:) by their corresponding global DOF indices 
% (assemblyPattern.u.R) and applies the @sum function. 
% This completely eliminates the need for slow, sequential for-loops 
% during the vector assembly.
R = accumarray( ...
    assemblyPattern.u.R, ...
    V_R(:), ...
    [ndof,1], ...
    @sum,0);

end