function bc = updateDisplacementBCquad(bc, lambda, load, mesh)
% ============================================================
% Updates displacement boundary values for the current
% normalized load parameter lambda.
%
% Physical displacement (Quadratic profile):
%       u_y(x) = lambda * maxDisplacement * f(x)
% ============================================================

% Initialize all constrained DOFs to zero
bc.u.constrainedValues = zeros(size(bc.u.constrainedDofs));

% Logical mask for constrained DOFs that are loaded
isLoaded = ismember(bc.u.constrainedDofs, bc.u.loadedDofs);

% Get the specific DOFs that are loaded
loaded_DOFs = bc.u.constrainedDofs(isLoaded);

% Assuming standard 2D DOF numbering (2 DOFs per node: u_x, u_y),
% the corresponding node index is ceil(DOF / 2).
nodeIndices = ceil(loaded_DOFs / 2);

% Extract x-coordinates for these nodes directly from the mesh structure
x = mesh.nodes.x(nodeIndices);

% Dynamically calculate plate width for normalization 
x_min = min(mesh.nodes.x);
x_max = max(mesh.nodes.x);
W = x_max - x_min;

% Normalize x to a [0, 1] range to easily apply the quadratic shape
x_norm = (x - x_min) / W;

% --- CHOOSE YOUR QUADRATIC PROFILE ---
% Option A: Max displacement at center, zero at edges (Parabolic)
quadraticProfile = 4 * x_norm .* (1 - x_norm); 

% Option B: Zero at left edge, max displacement at right edge (Asymmetric)
% quadraticProfile = x_norm.^2;

% Option C: Max at edges, zero at center (Inverted Parabola)
% quadraticProfile = 4 * (x_norm - 0.5).^2;

% Assign the prescribed non-uniform displacement 
bc.u.constrainedValues(isLoaded) = lambda * load.maxDisplacement .* quadraticProfile;

end