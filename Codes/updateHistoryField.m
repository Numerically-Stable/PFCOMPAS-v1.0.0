function Hnew = updateHistoryField( ...
    mesh,u,H_old,material,model,quad,parallel)
% ============================================================
% updateHistoryField - Enforce Thermodynamic Irreversibility
%
% PHASE-FIELD PHYSICS RELEVANCE:
% The second law of thermodynamics mandates that brittle fracture is an 
% irreversible dissipative process; cracks cannot physically heal. 
% Mathematically, this imposes a strict inequality constraint on the 
% damage variable: d(phi)/dt >= 0.
%
% To enforce this without resorting to computationally expensive active-set 
% or penalty methods within the phase-field Newton solver, the framework 
% implements the strain history approach proposed by Miehe et al. 
%
% The true active tensile strain energy (psi_plus) is replaced by the 
% history field (H), defined as the maximum tensile energy achieved over 
% the entire deformation history of the solid:
%       H_new = max(H_old, psi_plus_current)
%
% Because H is monotonically increasing, the energetic driving force for 
% damage can never decrease, naturally preventing crack healing during 
% compressive unloading.
%
% GAUSS POINT STORAGE:
% Notice that H is stored and evaluated strictly at the integration points 
% (nGP). Projecting this history field to the nodes via an L2 projection 
% would introduce artificial numerical diffusion (smearing), which severely 
% degrades the sharpness of the crack path.
% ============================================================

nelem = mesh.elem.n;
ngp   = numel(quad.W);

% Pre-allocate the updated history field array [nelem x ngp]
Hnew = zeros(nelem,ngp);

%% ============================================================
% Element Evaluation Loop
% ============================================================
% Evaluating the local strain energy is an embarrassingly parallel task.
if parallel.useParpool
    parfor e = 1:nelem
        % Compute the instantaneous active tensile strain energy 
        % (psi_plus) for all Gauss points in the current element based 
        % on the current displacement field u.
        psiPlus = ...
            elementPositiveEnergy( ...
            mesh,e,u,material,model,quad);
            
        % Update the history field by taking the element-wise maximum
        Hnew(e,:) = ...
            max(H_old(e,:),psiPlus');
    end
else
    for e = 1:nelem
        psiPlus = ...
            elementPositiveEnergy( ...
            mesh,e,u,material,model,quad);
            
        Hnew(e,:) = ...
            max(H_old(e,:),psiPlus');
    end
end

end