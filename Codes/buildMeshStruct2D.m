function mesh = buildMeshStruct2D(node, connectivity, elemType)
% buildMeshStruct2D
% ------------------------------------------------------------
% Build a structured mesh database for coupled FEM solvers (Q4)
%
% This function transforms standard raw mesh data into a highly 
% structured database. It precomputes separate Degree of Freedom (DOF) 
% maps for the multi-field problem (displacement vs. phase-field) 
% and automatically identifies geometric boundaries for Dirichlet/Neumann 
% boundary condition application.
%
% INPUT:
%   node         : [nnode x dim] global nodal coordinates
%   connectivity : [nelem x nen] global element connectivity matrix
%   elemType     : element type string (e.g. 'Q4')
%
% OUTPUT:
%   mesh : struct with fields containing topology, DOFs, and boundaries.
%
% Author: Pranjal Saxena 
% ------------------------------------------------------------

%% ------------------------------------------------------------
% Meta
%% ------------------------------------------------------------
mesh = struct();
mesh.meta.dim      = size(node,2);
mesh.meta.elemType = elemType;

% Enforce 4-node quadrilateral restriction. Q4 elements are optimal 
% for this framework's specific quadrature rules and error estimators.
if strcmp(elemType,'Q4')
    nen = 4;
else
    error('buildMeshStruct: only Q4 supported for now');
end

%% ------------------------------------------------------------
% Nodes
%% ------------------------------------------------------------
mesh.nodes = struct();
mesh.nodes.coord = node;
mesh.nodes.n     = size(node,1);
mesh.nodes.x     = node(:,1);
mesh.nodes.y     = node(:,2);

%% ------------------------------------------------------------
% Elements and Degree of Freedom (DOF) Mapping
%% ------------------------------------------------------------
mesh.elem = struct();
mesh.elem.conn = connectivity;
mesh.elem.n    = size(connectivity,1);
mesh.elem.nen  = nen;
mesh.elem.type = elemType;

% -------------------------------------------------------------------
% CRITICAL: Precompute DOF maps for the alternate minimization solver.
% The phase-field problem is a coupled multi-field problem. We must 
% map the local element nodes to the global system of equations 
% independently for the vector displacement field and scalar damage field.
% -------------------------------------------------------------------

% ---- Displacement DOFs (2 DOFs per node: ux, uy) ----
mesh.elem.dofU = zeros(mesh.elem.n, 2*nen);
mesh.elem.dofU(:,1:2:end) = 2*connectivity - 1;  % Odd indices: x-displacement (ux)
mesh.elem.dofU(:,2:2:end) = 2*connectivity;      % Even indices: y-displacement (uy)

% ---- Phase-field DOFs (1 DOF per node: phi) ----
mesh.elem.dofP = connectivity; % Direct 1-to-1 mapping for the scalar field

%% ------------------------------------------------------------
% Edges (topology)
%% ------------------------------------------------------------
% Call external helper to generate unique global edges and element mapping
edgeConnectivity = buildEdgeConnectivity(connectivity);
mesh.edge = struct();
mesh.edge.conn = edgeConnectivity;
mesh.edge.n    = size(edgeConnectivity,1);

% Local edge node numbering convention for standard Q4 isoparametric elements
localEdgeNodes = [1 2; 2 3; 3 4; 4 1];

% Extract the global node pairs for every edge in the mesh
mesh.edge.nodes = zeros(mesh.edge.n,2);
for k = 1:mesh.edge.n
    eL   = edgeConnectivity(k,1);
    locL = edgeConnectivity(k,2);
    mesh.edge.nodes(k,:) = connectivity(eL, localEdgeNodes(locL,:));
end

% Identify exterior boundary edges (where the element sharing count is zero/null)
mesh.edge.isBoundary = (mesh.edge.conn(:,3) == 0);

%% ------------------------------------------------------------
% Boundary classification (geometric extraction)
%% ------------------------------------------------------------
% Automates the identification of standard Cartesian boundaries to simplify 
% the application of boundary conditions (e.g., tensile loading, rigid supports).
mesh.boundary = struct();
tol = 1e-8; % Geometric tolerance to account for floating-point inaccuracies
y = mesh.nodes.y;
x = mesh.nodes.x;
edgeNodes = mesh.edge.nodes;
isBnd     = mesh.edge.isBoundary;

% Determine domain bounding box
ymax = max(y); ymin = min(y);
xmin = min(x); xmax = max(x);

% --- Top boundary (y = ymax) ---
isTop = isBnd & ...
        abs(y(edgeNodes(:,1)) - ymax) < tol & ...
        abs(y(edgeNodes(:,2)) - ymax) < tol;
mesh.boundary.top.edges = find(isTop);
mesh.boundary.top.nodes = unique(edgeNodes(isTop,:));

% --- Bottom boundary (y = ymin) ---
isBottom = isBnd & ...
           abs(y(edgeNodes(:,1)) - ymin) < tol & ...
           abs(y(edgeNodes(:,2)) - ymin) < tol;
mesh.boundary.bottom.edges = find(isBottom);
mesh.boundary.bottom.nodes = unique(edgeNodes(isBottom,:));

% --- Left boundary (x = xmin) ---
isLeft = isBnd & ...
         abs(x(edgeNodes(:,1)) - xmin) < tol & ...
         abs(x(edgeNodes(:,2)) - xmin) < tol;
mesh.boundary.left.edges = find(isLeft);
mesh.boundary.left.nodes = unique(edgeNodes(isLeft,:));

% --- Right boundary (x = xmax) ---
isRight = isBnd & ...
          abs(x(edgeNodes(:,1)) - xmax) < tol & ...
          abs(x(edgeNodes(:,2)) - xmax) < tol;
mesh.boundary.right.edges = find(isRight);
mesh.boundary.right.nodes = unique(edgeNodes(isRight,:));

%% ------------------------------------------------------------
% Sanity checks
%% ------------------------------------------------------------
mesh.meta.nnode = mesh.nodes.n;
mesh.meta.nelem = mesh.elem.n;
mesh.meta.nedge = mesh.edge.n;
fprintf('Mesh built: %d nodes, %d elements, %d edges\n', ...
        mesh.meta.nnode, mesh.meta.nelem, mesh.meta.nedge);
end