function R = computeReactionForce( ...
    mesh,u,phi,model,material,quad,boundaryName,dir)
% ==============================================================
% computeReactionForce - Integrated Boundary Traction Evaluation
%
% PHASE-FIELD PHYSICS RELEVANCE:
% To plot the macroscopic load-displacement response of the fracturing 
% structure, we must evaluate the total reaction force acting on the 
% driven boundary. Furthermore, this value is strictly required to 
% calculate the external work done (\Delta W_ext) for the global 
% thermodynamic energy balance check.
%
% Mathematically, this is the contour integral of the traction vector 
% along the specified Dirichlet boundary (Gamma):
%       R_dir = integral_Gamma (sigma * n)_dir dGamma
%
% HPC RELEVANCE:
% Because this is a post-processing step, we do not need to assemble 
% global arrays. We simply loop over the specific 1D boundary edges, 
% evaluate the degraded stress at the 1D edge Gauss points, and accumulate 
% the scalar force.
% ==============================================================

%% --------------------------------------------------------------
% Basic Information and Initialization
% --------------------------------------------------------------
R = 0;
dim = model.dim;
nen   = mesh.elem.nen;

% Standard local edge node mapping for Q4 elements
localEdgeNodes = [1 2; 2 3; 3 4; 4 1];

% 1D Edge quadrature rules (generated previously by Gauss1D.m)
We = quad.edgeW;
Qe = quad.edgeQ;

% Retrieve the pre-classified boundary edges from the mesh struct
edgeList = mesh.boundary.(boundaryName).edges(:);

%% ==============================================================
% Loop over boundary edges
% ==============================================================
for k = edgeList'
    % Extract global element ID and local edge ID
    eL   = mesh.edge.conn(k,1);
    locL = mesh.edge.conn(k,2);
    
    conn = mesh.elem.conn(eL,:);
    pts  = mesh.nodes.coord(conn,:);
    
    phi_e = phi(conn);
    dofU  = mesh.elem.dofU(eL,:);
    u_e   = u(dofU);
    
    % Extract physical coordinates of the two edge nodes
    edgeNodes = conn(localEdgeNodes(locL,:));
    x1 = mesh.nodes.coord(edgeNodes(1),:);
    x2 = mesh.nodes.coord(edgeNodes(2),:);
    
    % The 1D edge Jacobian maps the reference domain [-1, 1] to the 
    % physical edge length L. J_edge = L / 2.
    edgeJac = norm(x2 - x1)/2;
    
    % Tangent vector along the edge (from node 1 to node 2)
    t = (x2 - x1)/norm(x2 - x1);
    
    % Initial mathematical normal vector (90 degree clockwise rotation of t)
    n = [t(2); -t(1)];
    
    %% -------------------------------------------------------
    % CRITICAL FIX: Enforce Boundary-Based Outward Normal
    % -------------------------------------------------------
    % Depending on the global nodal numbering and connectivity winding 
    % (clockwise vs. counter-clockwise), the mathematically derived normal 
    % might point INTO the domain rather than OUTWARD. 
    % This block explicitly checks the assigned boundary orientation and 
    % flips the normal vector if necessary, strictly ensuring that the 
    % integrated traction has the correct thermodynamic sign.
    switch boundaryName
        case 'top'
            if n(2) < 0, n = -n; end
        case 'bottom'
            if n(2) > 0, n = -n; end
        case 'right'
            if n(1) < 0, n = -n; end
        case 'left'
            if n(1) > 0, n = -n; end
    end
    
    %% ==============================================================
    % 1D Edge Gauss Integration
    % ==============================================================
    for q = 1:length(We)
        s = Qe(q); % 1D coordinate in [-1, 1]
        w = We(q); % 1D integration weight
        
        % Map the 1D edge coordinate (s) to the 2D parent element 
        % coordinates (xi, eta).
        [xi, eta] = refEdgeCoords(locL, s);
        
        % Evaluate standard 2D shape functions at this edge point
        [N, dNdxi] = lagrange_basis(mesh.meta.elemType, [xi, eta]);
        
        J    = pts' * dNdxi;
        dNdx = dNdxi / J;
        
        %% --------------------------------------------------------
        % Kinematics
        % ---------------------------------------------------------
        B = displacementBMatrix(dNdx,dim,nen);
        strain = B * u_e;
        phi_gp = N' * phi_e;
        
        %% --------------------------------------------------------
        % Degradation and Energy Split
        % ---------------------------------------------------------
        % Even at the boundary, we must evaluate the exact degraded stress.
        [g,~,~] = degradationFunction(phi_gp,model,material);        
        split   = computeSplitInfo(strain,material,model);
    
        %% --------------------------------------------------------
        % Reconstruct Degraded Cauchy Stress Tensor
        % ---------------------------------------------------------
        switch lower(model.split)
            case 'none'
                sigma = g*split.sigma;
            case {'amor','spectral'}
                sigma = g*split.sigmaPlus + split.sigmaMinus;
            otherwise
                error('Unknown split: %s', model.split);
        end
    
        %% --------------------------------------------------------
        % Compute Traction Vector and Accumulate Force
        % ---------------------------------------------------------
        % Map 2nd-order stress tensor to 1D Voigt
        sigmaVoigt = stressTensorToVoigt(sigma,dim);
        
        % Reconstruct the 2x2 Cauchy stress matrix for matrix multiplication
        sigmaMat = [sigmaVoigt(1) sigmaVoigt(3);
                    sigmaVoigt(3) sigmaVoigt(2)];
                    
        % Cauchy's law: traction = sigma * n
        traction = sigmaMat * n;
        
        % Integrate: R = R + (traction_dir) * weight * Jacobian
        R = R + traction(dir) * w * edgeJac;
    end
end
end