function [Eelastic,Efracture,Etotal] = ...
    postprocessAT2Energy( ...
    mesh,u,phi,material,model,quad,parallel)
% ============================================================
% postprocessAT2Energy - Global Energy Integration
%
% PHASE-FIELD PHYSICS RELEVANCE:
% The core premise of the Francfort-Marigo variational fracture approach 
% is the minimization of a global energy functional. This functional 
% consists of two competing terms:
%   1. Eelastic: The global degraded elastic strain energy.
%   2. Efracture: The global dissipated fracture surface energy.
%
% HPC RELEVANCE:
% Calculating the exact global integrals at every single Newton-Raphson 
% iteration would impose a massive, unnecessary computational penalty. 
% Therefore, this postprocessing routine is strictly called ONLY after 
% the alternate minimization (staggered) solver has fully converged for 
% the current load step. 
%
% The resulting scalar energy values are critical for the global 
% thermodynamic balance check (evaluating whether External Work Done 
% equals the change in Internal Energy).
% ============================================================

nelem = mesh.elem.n;

% Initialize global scalar energy accumulators
Eelastic  = 0;
Efracture = 0;

%% ============================================================
% Element Evaluation Loop (Map-Reduce Architecture)
% ============================================================
% Because the global energy is simply the spatial integral of the local 
% energy densities, the element-wise evaluations are completely independent 
% and highly parallelizable.
if parallel.useParpool
    % Pre-allocate arrays to hold the scalar energy of each element
    EelElem = zeros(nelem,1);
    EfElem  = zeros(nelem,1);
    
    parfor e = 1:nelem
        % Compute the integrated energy for the current element
        [EelElem(e),EfElem(e)] = ...
            elementEnergy( ...
            mesh,e,u,phi, ...
            material,model,quad);
    end
    
    % Map-Reduce: Sum the element arrays into global scalars
    Eelastic  = sum(EelElem);
    Efracture = sum(EfElem);
else
    for e = 1:nelem
        [Eel,Ef] = ...
            elementEnergy( ...
            mesh,e,u,phi, ...
            material,model,quad);
            
        Eelastic  = Eelastic  + Eel;
        Efracture = Efracture + Ef;
    end
end

%% ============================================================
% Total Potential Energy
% ============================================================
Etotal = Eelastic + Efracture;

end