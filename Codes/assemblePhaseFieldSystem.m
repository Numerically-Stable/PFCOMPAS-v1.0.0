function [Kphi,Rphi] = ...
    assemblePhaseFieldSystem( ...
    mesh,phi,H,model,material,quad, ...
    parallel,assemblyPattern)
% ============================================================
% assemblePhaseFieldSystem - Global Phase-Field Assembly
%
% PHASE-FIELD HPC RELEVANCE:
% Much like the displacement Newton solver, the phase-field Newton-Raphson 
% loop requires the continuous reassembly of the global tangent stiffness 
% matrix (Kphi) and the internal residual force vector (Rphi). 
%
% Because the phase-field equation is heavily non-linear with respect to 
% phi (due to the degradation derivatives g'(phi) and g''(phi)), this 
% assembly is the primary bottleneck for the damage evolution step.
%
% This subroutine strictly replicates the high-performance map-reduce 
% paradigm established in assembleDisplacementSystem. It decouples the 
% mathematical integration from the sparse allocation, utilizing PARFOR 
% for thread-safe concurrent element evaluation, followed by vectorized 
% sparse() and accumarray() calls for the global matrix construction.
%
% Note: The phase-field is a scalar variable (1 DOF per node), making 
% these assembly arrays significantly smaller than their mechanical 
% counterparts.
% ============================================================

nelem  = mesh.elem.n;
nen    = mesh.elem.nen;
nnodes = mesh.nodes.n;

%% ============================================================
% Element numerical values (Pre-allocation)
% ============================================================
% V_K holds the flattened local tangent stiffness matrices.
% V_R holds the local residual force vectors.
% Size is based strictly on the scalar nodes per element (nen).
V_K = zeros(nen*nen,nelem);
V_R = zeros(nen,nelem);

%% ============================================================
% Element loop (Integration Phase)
% ============================================================
if parallel.useParfor
    parfor e = 1:nelem
        % Compute local tangent stiffness (Ke) and internal forces (Re)
        % based on the current damage state (phi) and history field (H).
        [Ke,Re] = ...
            elementPhaseFieldContribution( ...
            mesh,e,phi,H,model,material,quad);
            
        % Store values in independent column slices to guarantee thread safety
        V_K(:,e) = Ke(:);
        V_R(:,e) = Re;
    end
else
    for e = 1:nelem
        [Ke,Re] = ...
            elementPhaseFieldContribution( ...
            mesh,e,phi,H,model,material,quad);
        V_K(:,e) = Ke(:);
        V_R(:,e) = Re;
    end
end

%% ============================================================
% Sparse global stiffness (Tangent Matrix Kphi)
% ============================================================
% Vectorized summation of overlapping node stiffness contributions.
Kphi = sparse( ...
    assemblyPattern.phi.I, ...
    assemblyPattern.phi.J, ...
    V_K(:), ...
    nnodes,nnodes);

%% ============================================================
% Global residual (Internal Force Vector Rphi)
% ============================================================
% Vectorized assembly of the 1D global residual vector, mathematically 
% identical to the mechanical assembly logic.
Rphi = accumarray( ...
    assemblyPattern.phi.R, ...
    V_R(:), ...
    [nnodes,1], ...
    @sum,0);

end