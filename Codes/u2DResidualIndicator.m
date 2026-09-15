function [Eta,Eta_elem,Eta_elem2] = ...
    u2DResidualIndicator( ...
    mesh,u,phi,material,model,quad,dir_nodes)
% ============================================================
% u2DResidualIndicator - Verfürth A Posteriori Error Estimator (Mechanics)
%
% PHASE-FIELD HPC & MATHEMATICAL RELEVANCE:
% Standard Newton-Raphson residuals only tell us if the algebraic equations 
% (K*du = -R) have been solved correctly. They tell us NOTHING about the 
% discretization error (i.e., the difference between the exact continuous 
% solution and our discrete finite element approximation).
%
% To gauge the spatial accuracy of the mesh and the length-scale l0, we 
% employ explicit a posteriori error estimators based on Verfürth's formulations.
% For the mechanical equilibrium equation div(g(phi)*sigma) = 0, the error 
% is bounded by two terms:
%   1. Element Volumetric Residual: Forces violating equilibrium inside the element.
%   2. Edge Jump Residual: Discontinuities in the traction vector across 
%      element boundaries, [t] = (sigma_L - sigma_R)*n.
%
% NOTE ON DIRICHLET BOUNDARIES:
% Error estimators evaluate the violation of the weak form. Because displacements 
% are exactly prescribed on Dirichlet boundaries, the test functions delta_u 
% are zero there. Thus, Dirichlet edges are skipped.
% ============================================================

nelem = mesh.elem.n;
nen   = mesh.elem.nen;
dim   = model.dim;

% Initialize squared element indicator array
Eta_elem2 = zeros(nelem,1);

% Volumetric Quadrature
W = quad.W;
Q = quad.Q;

%% ============================================================
% (1) ELEMENT VOLUMETRIC RESIDUALS
%% ============================================================
for e = 1:nelem
    conn  = mesh.elem.conn(e,:);
    pts   = mesh.nodes.coord(conn,:);
    
    volAccum = 0;
    for q = 1:length(W)
        pt = Q(q,:);
        wt = W(q);
        [~, dNdxi] = lagrange_basis(mesh.meta.elemType, pt);
        J    = pts' * dNdxi;
        detJ = det(J);
        
        % ----------------------------------------------------
        % Q4 Divergence Approximation
        % ----------------------------------------------------
        % For standard low-order Q4 bilinear elements, the second spatial 
        % derivatives of the shape functions evaluate to zero (or are 
        % vanishingly small). Therefore, div(g*sigma) inside the element 
        % is mathematically negligible compared to the edge flux jumps. 
        % Without advanced gradient reconstruction (e.g., SPR), this 
        % term is approximated as zero.
        divSigma = [0;0];
        
        volAccum = volAccum + (norm(divSigma)^2) * wt * detJ;
    end
    
    % Element characteristic size (hK)
    area = polyarea(pts(:,1), pts(:,2));
    hK   = sqrt(area);
    
    % Scale the volumetric residual by the square of the element size
    Eta_elem2(e) = (hK^2) * volAccum;
end

%% ============================================================
% (2) EDGE JUMP TERMS (TRACTION DISCONTINUITY)
%% ============================================================
% The dominant source of discretization error in low-order FEM.
We = quad.edgeW;
Qe = quad.edgeQ;
localEdgeNodes = [1 2; 2 3; 3 4; 4 1];

for k = 1:mesh.edge.n
    edgeNodes = mesh.edge.nodes(k,:);
    
    % --------------------------------------------------------
    % Dirichlet Masking
    % --------------------------------------------------------
    % Skip edges entirely constrained by Dirichlet boundary conditions.
    if all(ismember(edgeNodes, dir_nodes))
        continue
    end
    
    % Extract Left and Right elements sharing this edge
    eL   = mesh.edge.conn(k,1);
    locL = mesh.edge.conn(k,2);
    eR   = mesh.edge.conn(k,3);
    
    connL = mesh.elem.conn(eL,:);
    ptsL  = mesh.nodes.coord(connL,:);
    dofUL = mesh.elem.dofU(eL,:);
    uL    = u(dofUL);
    phiL  = phi(connL);
    
    hasRight = (eR > 0);
    if hasRight
        connR = mesh.elem.conn(eR,:);
        ptsR  = mesh.nodes.coord(connR,:);
        dofUR = mesh.elem.dofU(eR,:);
        uR    = u(dofUR);
        phiR  = phi(connR);
    end
    
    % --------------------------------------------------------
    % Edge Geometry and Normal
    % --------------------------------------------------------
    n1 = localEdgeNodes(locL,1);
    n2 = localEdgeNodes(locL,2);
    x1 = ptsL(n1,:);
    x2 = ptsL(n2,:);
    
    edgeLength = norm(x2 - x1);
    edgeJac = edgeLength/2;
    t = (x2 - x1)/edgeLength;
    n = [t(2); -t(1)]; % Normal vector pointing from Left to Right
    
    %% --- 1D Edge Gauss Integration ---
    for q = 1:length(We)
        s = Qe(q);
        w = We(q);
        
        %% -------- EVALUATE LEFT ELEMENT TRACTION --------
        [xiL, etaL] = refEdgeCoords(locL,s);
        [N_L, dNdxiL] = lagrange_basis(mesh.meta.elemType,[xiL,etaL]);
        J0L = ptsL' * dNdxiL;
        dNdxL = dNdxiL / J0L;
    
        B_L = displacementBMatrix(dNdxL,dim,nen);
        strainL = B_L * uL;
        phiL_q = N_L' * phiL;
    
        [g,~,~] = degradationFunction(phiL_q,model,material);
        split   = computeSplitInfo(strainL,material,model);
    
        switch lower(model.split)
            case 'none'
                sigma = g*split.sigma;
            case {'amor','spectral'}
                sigma = g*split.sigmaPlus + split.sigmaMinus;
            otherwise
                error('Unknown split: %s', model.split);
        end
    
        sigmaL = stressTensorToVoigt(sigma,dim);
        tractionL = [ sigmaL(1), sigmaL(3);
                      sigmaL(3), sigmaL(2) ] * n;
                           
        %% -------- EVALUATE RIGHT ELEMENT TRACTION --------
        if hasRight
            locR = mesh.edge.conn(k,4);
            [xiR, etaR] = refEdgeCoords(locR,s);
            [N_R, dNdxiR] = lagrange_basis(mesh.meta.elemType,[xiR,etaR]);
            J0R = ptsR' * dNdxiR;
            dNdxR = dNdxiR / J0R;
        
            B_R = displacementBMatrix(dNdxR,dim,nen);
            strainR = B_R * uR;           
            phiR_q = N_R' * phiR;
        
            [g,~,~] = degradationFunction(phiR_q,model,material);
            split = computeSplitInfo(strainR,material,model);
        
            switch lower(model.split)
                case 'none'
                    sigma = g*split.sigma;
                case {'amor','spectral'}
                    sigma = g*split.sigmaPlus + split.sigmaMinus;
                otherwise
                    error('Unknown split: %s', model.split);
            end
        
            sigmaR = stressTensorToVoigt(sigma,dim);
            tractionR =  [ sigmaR(1), sigmaR(3);
                           sigmaR(3), sigmaR(2) ] * n;
        else
            % For an unconstrained Neumann boundary, traction should be 0.
            tractionR = [0;0];
        end
        
        %% -------- EVALUATE FLUX JUMP --------
        % [t] = t_L - t_R
        jump = tractionL - tractionR;
        
        % The edge contribution is scaled by the element edge length (h_e)
        contrib = edgeLength * dot(jump,jump) * w * edgeJac;
        
        % Distribute the error equally to both adjacent elements
        if hasRight
            Eta_elem2(eL) = Eta_elem2(eL) + 0.5*contrib;
            Eta_elem2(eR) = Eta_elem2(eR) + 0.5*contrib;
        else
            Eta_elem2(eL) = Eta_elem2(eL) + contrib;
        end
    end
end

%% ============================================================
% FINALIZE INDICATORS
%% ============================================================
% --- Detect NaNs in squared indicator ---
nan_idx = find(isnan(Eta_elem2));
if ~isempty(nan_idx)
    warning('NaNs detected in Eta_elem2 at elements: %s', mat2str(nan_idx));
end

% --- Clean NaNs (prevent propagation) ---
Eta_elem2_clean = Eta_elem2;
Eta_elem2_clean(isnan(Eta_elem2_clean)) = 0;

% --- Sanity check ---
if any(Eta_elem2_clean < 0)
    error('Negative values detected in Eta_elem2 (should not happen)');
end

% --- Element-wise indicator ---
Eta_elem = sqrt(Eta_elem2_clean);

% --- Global indicator ---
Eta = sqrt(sum(Eta_elem2_clean));

end