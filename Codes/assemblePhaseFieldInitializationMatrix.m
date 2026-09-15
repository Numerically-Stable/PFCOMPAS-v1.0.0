function A = ...
    assemblePhaseFieldInitializationMatrix( ...
    mesh,material,quad,parallel,assemblyPattern)
% ============================================================
% assemblePhaseFieldInitializationMatrix
%
% This function constructs the global initialization matrix:
%     A = M + l0^2 K
% corresponding to the homogeneous screened Poisson equation:
%     phi - l0^2 Laplacian(phi) = 0
%
% PHASE-FIELD HPC RELEVANCE:
% To achieve maximum performance, this routine strictly decouples the 
% element-level numerical integration from the global sparse matrix assembly. 
% By pre-allocating a dense array (V_K) to catch the flattened local matrices, 
% we can utilize MATLAB's parallel computing toolbox (parfor) to evaluate 
% elements concurrently without triggering memory race conditions or 
% requiring complex atomic locks.
% ============================================================

nelem  = mesh.elem.n;
nen    = mesh.elem.nen;
nnodes = mesh.nodes.n;

%% ============================================================
% Element numerical values pre-allocation
% ============================================================
% V_K holds the flattened local stiffness matrices for every element.
% By sizing it as [nen*nen, nelem], each column corresponds to exactly 
% one element. This specific memory layout is required for safe parfor execution.
V_K = zeros(nen*nen,nelem);

%% ============================================================
% Element loop (Integration Phase)
% ============================================================
% The parallel pool distributes chunks of elements to different workers.
if parallel.useParfor
    parfor e = 1:nelem
        % Compute the local initialization matrix (M_e + l0^2 * K_e)
        Ke = elementPhaseFieldInitialization( ...
            mesh, ...
            e, ...
            material, ...
            quad);
            
        % Store the flattened matrix into the independent column slice.
        % This completely avoids race conditions between parallel threads.
        V_K(:,e) = Ke(:);
    end
else
    % Fallback to serial execution if parpool is unavailable or disabled.
    for e = 1:nelem
        Ke = elementPhaseFieldInitialization( ...
            mesh, ...
            e, ...
            material, ...
            quad);
        V_K(:,e) = Ke(:);
    end
end

%% ============================================================
% Global sparse matrix (Assembly Phase)
% ============================================================
% MATLAB's sparse() function natively handles the summation of overlapping 
% DOFs. If multiple elements contribute to the same (I, J) node pair, 
% sparse() automatically accumulates (adds) the values from V_K(:). 
% This single vectorized call replaces the traditional, painfully slow 
% nested assembly loops.
A = sparse( ...
    assemblyPattern.phi.I, ...
    assemblyPattern.phi.J, ...
    V_K(:), ...
    nnodes, ...
    nnodes);
    
%% ============================================================
% Enforce strict numerical symmetry
% ============================================================
% While finite element matrices are analytically symmetric, independent 
% floating-point operations across different parallel workers can introduce 
% vanishingly small asymmetries (on the order of 10^-16). 
% If 'A' is not strictly symmetric in memory, efficient direct solvers 
% (like Cholesky decomposition) will fail or silently fall back to slower, 
% generic LU decompositions. This explicit symmetrization prevents that.
A = 0.5*(A + A.');

end