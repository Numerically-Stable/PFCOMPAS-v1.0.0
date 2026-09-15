function Hnew = ...
    calculateHistory( ...
    mesh,u,Hold,material,model,quad,parallel)
% ============================================================
% calculateHistory
%
% Updates the history field:
%
%       Hnew = max(Hold, psiPlus)
%
% This function is called after the displacement Newton solve
% has converged.
%
% Only Hnew is retained globally.
% No GP stresses or tangent matrices are stored.
%
% ============================================================

nelem = mesh.elem.n;

ngp = numel(quad.W);


%% ============================================================
% Allocate history
% ============================================================

Hnew = zeros(nelem,ngp);


%% ============================================================
% Element loop
% ============================================================

if parallel.useParfor

    parfor e = 1:nelem

        psiPlus = ...
            elementPositiveEnergy( ...
            mesh,e,u,material,model,quad);

        Hnew(e,:) = ...
            max(Hold(e,:),psiPlus(:)');

    end

else

    for e = 1:nelem

        psiPlus = ...
            elementPositiveEnergy( ...
            mesh,e,u,material,model,quad);

        Hnew(e,:) = ...
            max(Hold(e,:),psiPlus(:)');

    end

end

end