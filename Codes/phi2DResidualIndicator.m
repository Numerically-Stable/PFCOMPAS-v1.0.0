function [Eta,Eta_elem,Eta_elem2] = ...
    phi2DResidualIndicator( ...
    mesh,phi,H,material,model,quad,crack_nodes)
% ------------------------------------------------------------
% phi2DResidualIndicator - Verfürth A Posteriori Error Estimator (Damage)
%
% PHASE-FIELD HPC & MATHEMATICAL RELEVANCE:
% Similar to the mechanical equilibrium problem, the spatial discretization 
% of the phase-field PDE introduces inherent errors. The exact continuous 
% strong form of the AT2 damage evolution equation is:
%
%       g'(phi)*H + (Gc/l0)*phi - Gc*l0*Laplacian(phi) = 0
%
% To estimate the discretization error (Eta_phi), we must measure the 
% continuous violation of this PDE over the finite element mesh. This is 
% bounded by:
%   1. The volumetric residual (forces attempting to evolve damage internally).
%   2. The edge flux jump (discontinuities in the phase-field gradient 
%      across element boundaries).
%
% DIRICHLET MASKING:
% The nodes defining the initial macroscopic discrete crack (crack_nodes) 
% are heavily constrained (phi = 1.0). Edges lying entirely on this 
% boundary do not contribute to the weak-form test space and are 
% appropriately skipped.
% ------------------------------------------------------------

nelem = mesh.elem.n;
Eta_elem2 = zeros(nelem,1);

Gc = material.Gc;
l0 = material.l0;

W = quad.W;
Q = quad.Q;

%% ============================================================
% (1) ELEMENT VOLUMETRIC RESIDUALS
%% ============================================================
for e = 1:nelem
    conn  = mesh.elem.conn(e,:);
    pts   = mesh.nodes.coord(conn,:);
    phi_e = phi(conn);
    
    % The history field H is retrieved strictly at the Gauss points 
    % to prevent numerical diffusion/smearing.
    H_e   = H(e,:);   
    
    volAccum = 0;
    for q = 1:length(W)
        pt = Q(q,:);
        wt = W(q);
        
        [N, dNdxi] = lagrange_basis(mesh.meta.elemType, pt);
        J  = pts' * dNdxi;
        detJ = det(J);
        
        phi_gp = N' * phi_e;
    
        %% --------------------------------------------------------
        % Degradation Driving Force
        % ---------------------------------------------------------
        [~,dg,~] = degradationFunction(phi_gp,model,material);
        
        %% --------------------------------------------------------
        % Evaluate the Strong Form PDE
        % ---------------------------------------------------------
        % For standard bilinear Q4 elements, the Laplacian of phi evaluates 
        % to zero inside the element interior. The remaining terms are the 
        % thermodynamic driving force (dg * H) and the geometric restoration 
        % force (Gc/l0 * phi).
        Rphi_q = dg * H_e(q) + (Gc/l0)*phi_gp;
        
        volAccum = volAccum + (Rphi_q^2) * wt * detJ;
    end
    
    % Element characteristic size (area scaling)
    area = polyarea(pts(:,1), pts(:,2));
    hK = sqrt(area);
    
    % Volumetric error contribution scaled by the element size squared
    Eta_elem2(e) = (hK^2) * volAccum;
end

%% ============================================================
% (2) EDGE JUMP TERMS (GRADIENT FLUX DISCONTINUITY)
%% ============================================================
% The discontinuous gradient flux across the Q4 boundaries constitutes 
% the majority of the phase-field discretization error.
We = quad.edgeW;
Qe = quad.edgeQ;
localEdgeNodes = [1 2; 2 3; 3 4; 4 1];

for k = 1:mesh.edge.n
    edgeNodes = mesh.edge.nodes(k,:);
    
    % --------------------------------------------------------
    % Dirichlet Masking (Initial Cracks)
    % --------------------------------------------------------
    % Skip internal crack edges where phi is explicitly constrained.
    if all(ismember(edgeNodes, crack_nodes))
        continue
    end
    
    eL   = mesh.edge.conn(k,1);
    locL = mesh.edge.conn(k,2);
    eR   = mesh.edge.conn(k,3);
    
    connL = mesh.elem.conn(eL,:);
    ptsL  = mesh.nodes.coord(connL,:);
    phiL  = phi(connL);
    
    hasRight = (eR > 0);
    if hasRight
        connR = mesh.elem.conn(eR,:);
        ptsR  = mesh.nodes.coord(connR,:);
        phiR  = phi(connR);
    end
    
    % Edge geometry and outward normal vector (Left to Right)
    n1 = localEdgeNodes(locL,1);
    n2 = localEdgeNodes(locL,2);
    x1 = ptsL(n1,:);
    x2 = ptsL(n2,:);
    
    edgeLength = norm(x2 - x1);
    t = (x2 - x1)/edgeLength;
    n = [t(2); -t(1)];
    
    %% --- 1D Edge Gauss Integration ---
    for q = 1:length(We)
        s = Qe(q);
        w = We(q);
        
        %% -------- EVALUATE LEFT ELEMENT FLUX --------
        [xiL, etaL] = refEdgeCoords(locL, s);
        [~, dNdxiL] = lagrange_basis(mesh.meta.elemType,[xiL,etaL]);
        J0L   = ptsL' * dNdxiL;
        dNdxL = dNdxiL / J0L;
        
        gradphiL = dNdxL' * phiL;
        
        % The phase-field flux vector is F = -Gc * l0 * grad(phi).
        % Projected onto the normal: F * n
        fluxL = -Gc*l0 * (gradphiL' * n);
        
        %% -------- EVALUATE RIGHT ELEMENT FLUX --------
        if hasRight
            locR = mesh.edge.conn(k,4);
            [xiR, etaR] = refEdgeCoords(locR, s);
            [~, dNdxiR] = lagrange_basis(mesh.meta.elemType,[xiR,etaR]);
            J0R   = ptsR' * dNdxiR;
            dNdxR = dNdxiR / J0R;
            
            gradphiR = dNdxR' * phiR;
            fluxR = -Gc*l0 * (gradphiR' * n);
        else
            fluxR = 0; % Unconstrained Neumann boundary condition
        end
        
        %% -------- EVALUATE FLUX JUMP --------
        % [F] = F_L - F_R
        jump = fluxL - fluxR;
        
        % Verfürth scaling: the edge contribution is scaled by the 
        % physical length of the edge.
        contrib = edgeLength * (jump^2) * w * (edgeLength/2);
        
        % Distribute the error penalty symmetrically to the adjacent elements
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
% --- Detect NaNs in raw squared indicator ---
nan_idx = find(isnan(Eta_elem2));
if ~isempty(nan_idx)
    warning('NaNs detected in Eta_elem2 at elements: %s', mat2str(nan_idx));
end

% --- Clean NaNs (do NOT propagate into norms) ---
Eta_elem2_clean = Eta_elem2;
Eta_elem2_clean(isnan(Eta_elem2_clean)) = 0;

% --- Element-wise phase-field indicator ---
Eta_elem = sqrt(Eta_elem2_clean);

% --- Global phase-field indicator ---
Eta = sqrt(sum(Eta_elem2_clean));

end