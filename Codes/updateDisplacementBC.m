function bc = updateDisplacementBC(bc,lambda,load)
% ============================================================
% updateDisplacementBC - Incremental Dirichlet Boundary Conditions
%
% COMPUTATIONAL MECHANICS RELEVANCE:
% Phase-field fracture is an irreversible, path-dependent dissipative 
% process. To accurately track the evolution of the crack and capture 
% the highly non-linear softening behavior of the structure, the external 
% mechanical loading must be applied incrementally.
%
% This framework utilizes a displacement-controlled loading scheme, governed 
% by a monotonically increasing pseudo-time parameter, lambda, where:
%       lambda = 0.0 (Unloaded, pristine state)
%       lambda = 1.0 (Maximum prescribed macroscopic displacement)
%
% At each adaptive load step, this subroutine evaluates the current 
% lambda value and precisely updates the boundary condition array. It 
% carefully distinguishes between actively driven DOFs (the loading surface) 
% and statically fixed DOFs (the structural supports).
% ============================================================

% ------------------------------------------------------------
% Initialize all constrained DOFs to zero
% ------------------------------------------------------------
% This safely resets the constraint array, ensuring that rigid supports 
% (e.g., pinned bases or symmetry planes) are strictly locked to a 
% displacement of 0.0 at every load step.
bc.u.constrainedValues = zeros(size(bc.u.constrainedDofs));

% ------------------------------------------------------------
% Isolate and update the actively loaded DOFs
% ------------------------------------------------------------
% We use ismember to find the exact indices within the constrained array 
% that correspond to the actively driven nodes (e.g., the top edge of a 
% SENT specimen).
isLoaded = ismember(bc.u.constrainedDofs, bc.u.loadedDofs);

% ------------------------------------------------------------
% Apply the incremental displacement
% ------------------------------------------------------------
% The prescribed physical displacement vector is linearly scaled by lambda.
% This mapped value is strictly enforced during the Newton-Raphson assembly
% to drive the structural deformation and the subsequent strain history.
bc.u.constrainedValues(isLoaded) = lambda * load.maxDisplacement;

end