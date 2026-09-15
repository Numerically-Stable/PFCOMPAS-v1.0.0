function phi = initializePhaseField( ...
    mesh,crackNodes,material,quad,parallel,assemblyPattern)
% ============================================================
% initializePhaseField - Compute initial diffuse crack profile
%
% PHASE-FIELD PHYSICS RELEVANCE:
% Standard finite element meshes often define initial cracks as discrete, 
% sharp entities (e.g., a line of nodes where damage is assumed to be 1). 
% However, the AT2 phase-field model requires a continuous, regularized 
% damage band of width proportional to the length scale (l0). 
%
% To physically transition the sharp discrete crack into a mathematically 
% valid initial phase-field profile, we must solve the homogeneous 
% screened Poisson equation (which is the Euler-Lagrange equation of the 
% fracture surface functional in the absence of elastic strain energy):
%
%       phi - l0^2 * Laplacian(phi) = 0
%
% subject to the strict Dirichlet boundary condition:
%
%       phi = 1   on the initial discrete crack nodes.
%
% The resulting discrete finite element system is:
%
%       A * phi = 0
%
% where A = M + (l0^2)*K. 
% (M is the mass matrix and K is the phase-field diffusion/Laplacian matrix).
% ============================================================

ndof = mesh.nodes.n;

%% ============================================================
% Validate crack nodes
% ============================================================
crackNodes = unique(crackNodes(:));
if isempty(crackNodes)
    error('initializePhaseField:NoCrackNodes', ...
          'No initial crack nodes were supplied.');
end
if any(crackNodes < 1) || any(crackNodes > ndof)
    error('initializePhaseField:InvalidCrackNodes', ...
          'Invalid initial crack node index.');
end

%% ============================================================
% Assemble initialization matrix (A = M + l0^2 * K)
% ============================================================
% Calls a dedicated sub-routine to build the linear, load-independent 
% matrix 'A' without any elastic coupling or degradation functions.
A = assemblePhaseFieldInitializationMatrix( ...
    mesh, material, quad, parallel, assemblyPattern);

%% ============================================================
% Apply Dirichlet Boundary Conditions (Partitioning)
% ============================================================
fixedDofs = crackNodes;
isFixed = false(ndof,1);
isFixed(fixedDofs) = true;
freeDofs = find(~isFixed);

%% ============================================================
% Initialize phase field array
% ============================================================
phi = zeros(ndof,1);
phi(fixedDofs) = 1.0; % Enforce fully broken state exactly at crack nodes

%% ============================================================
% Solve the reduced linear algebraic system
%
% By partitioning the global matrix A into free (f) and constrained (c) 
% degrees of freedom, the homogeneous system (A * phi = 0) becomes:
%
%       [ Aff  Afc ] { phi_f } = { 0 }
%       [ Acf  Acc ] { phi_c }   { 0 }
%
% Which yields the reduced solvable system:
%       Aff * phi_f = -Afc * phi_c
% ============================================================
if ~isempty(freeDofs)
    Aff = A(freeDofs, freeDofs);
    Afc = A(freeDofs, fixedDofs);
    
    % The right-hand side is purely driven by the prescribed phi=1 boundary
    rhs = -Afc * phi(fixedDofs);
    
    % Direct linear solve for the unknown free DOFs
    phi(freeDofs) = Aff \ rhs;
end

%% ============================================================
% Enforce crack values exactly (Sanity Check)
% ============================================================
phi(fixedDofs) = 1.0;

%% ============================================================
% Numerical bounds and Floating-Point Correction
% ============================================================
% The continuous solution is strictly bounded between 0 and 1. 
% Small finite element interpolation undershoots or overshoots due to 
% machine precision are clamped here to prevent non-physical damage states.
phi(phi < 0 & phi > -1e-12) = 0;
phi(phi > 1 & phi < 1+1e-12) = 1;

end