function assemblyPattern = buildAssemblyPattern(mesh)
% ============================================================
% buildAssemblyPattern - Precompute sparse matrix allocation indices
%
% PHASE-FIELD HPC RELEVANCE:
% In the staggered phase-field solver, the global stiffness matrices 
% for both the displacement field (Ku) and the phase-field (Kphi) must 
% be reassembled hundreds of times per load step due to the highly 
% non-linear degradation function g(phi) and the Newton-Raphson iterations.
%
% Dynamically inserting local element matrices into a global sparse matrix 
% inside a loop (e.g., K(sctr, sctr) = K(sctr, sctr) + Ke) is severely 
% inefficient in MATLAB. To achieve near-C++ assembly speeds, we must use 
% the built-in sparse(I, J, V) function. 
% 
% This routine precomputes the global row (I) and column (J) index arrays 
% once during preprocessing. Because it strictly computes topological mapping, 
% no geometric or constitutive quantities (B-matrices, Jacobians) are stored 
% here, minimizing memory footprint.
%
% INPUT:
%   mesh : Unified mesh database struct
%
% OUTPUT:
%   assemblyPattern : Struct containing the I, J indices for Ku and Kphi, 
%                     and R indices for the residual vectors.
% ============================================================

nelem  = mesh.elem.n;
nnodes = mesh.nodes.n;
nen    = mesh.elem.nen;
ndofu  = size(mesh.elem.dofU, 2);

%% ============================================================
% DISPLACEMENT PATTERN (u-field)
% ============================================================
% Preallocate arrays based on exactly how many elements will contribute 
% to the global displacement stiffness matrix and residual vector.
nKu = nelem * ndofu * ndofu;
nRu = nelem * ndofu;

assemblyPattern.u.I  = zeros(nKu, 1);
assemblyPattern.u.J  = zeros(nKu, 1);
assemblyPattern.u.R  = zeros(nRu, 1);

for e = 1:nelem
    sctr = mesh.elem.dofU(e, :);
    
    %% --------------------------------------------------------
    % Stiffness Indices (I, J)
    %
    % ndgrid rapidly creates the Cartesian product of the local DOF 
    % scatter vector. 
    % I(row,col) = global row DOF
    % J(row,col) = global column DOF
    % This ordering perfectly matches the column-major flattening 
    % of the local elemental stiffness matrix Ke(:).
    % ---------------------------------------------------------
    [ii, jj] = ndgrid(sctr, sctr);
    
    kStart = (e - 1) * ndofu * ndofu + 1;
    kEnd   = e * ndofu * ndofu;
    
    assemblyPattern.u.I(kStart:kEnd) = ii(:);
    assemblyPattern.u.J(kStart:kEnd) = jj(:);
    
    %% --------------------------------------------------------
    % Residual Indices (R)
    % ---------------------------------------------------------
    rStart = (e - 1) * ndofu + 1;
    rEnd   = e * ndofu;
    assemblyPattern.u.R(rStart:rEnd) = sctr(:);
end

%% ============================================================
% PHASE-FIELD PATTERN (phi-field)
% ============================================================
% The phase-field is a scalar variable, requiring only 1 DOF per node. 
% Therefore, the local matrices are smaller (nen x nen).
nKphi = nelem * nen * nen;
nRphi = nelem * nen;

assemblyPattern.phi.I = zeros(nKphi, 1);
assemblyPattern.phi.J = zeros(nKphi, 1);
assemblyPattern.phi.R = zeros(nRphi, 1);

for e = 1:nelem
    sctr = mesh.elem.conn(e, :); % Direct 1-to-1 nodal mapping for phi
    
    %% --------------------------------------------------------
    % Stiffness Indices (I, J)
    % ---------------------------------------------------------
    [ii, jj] = ndgrid(sctr, sctr);
    
    kStart = (e - 1) * nen * nen + 1;
    kEnd   = e * nen * nen;
    
    assemblyPattern.phi.I(kStart:kEnd) = ii(:);
    assemblyPattern.phi.J(kStart:kEnd) = jj(:);
    
    %% --------------------------------------------------------
    % Residual Indices (R)
    % ---------------------------------------------------------
    rStart = (e - 1) * nen + 1;
    rEnd   = e * nen;
    assemblyPattern.phi.R(rStart:rEnd) = sctr(:);
end

%% ============================================================
% Store dimensions for rapid retrieval during active solve
% ============================================================
assemblyPattern.nnode = nnodes;
assemblyPattern.nelem = nelem;
assemblyPattern.nen   = nen;
assemblyPattern.ndofu = ndofu;
end