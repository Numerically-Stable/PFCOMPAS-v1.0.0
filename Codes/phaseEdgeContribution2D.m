function [eL,eR,contribution] = ...
    phaseEdgeContribution2D( ...
    mesh,edgeID,phi,material,model,quad)
% ==============================================================
% PHASE-FIELD FLUX JUMP ON ONE EDGE
%
% Returns:
%
%   eL
%   eR
%   contribution
%
% where contribution is the complete Verfürth edge term.
% ==============================================================


%% --------------------------------------------------------------
% Edge connectivity
% --------------------------------------------------------------

eL   = mesh.edge.conn(edgeID,1);
locL = mesh.edge.conn(edgeID,2);

eR = mesh.edge.conn(edgeID,3);

hasRight = eR > 0;


%% --------------------------------------------------------------
% Left element
% --------------------------------------------------------------

connL = mesh.elem.conn(eL,:);

ptsL = mesh.nodes.coord(connL,:);

phiL = phi(connL);


%% --------------------------------------------------------------
% Right element
% --------------------------------------------------------------

if hasRight

    locR = mesh.edge.conn(edgeID,4);

    connR = mesh.elem.conn(eR,:);

    ptsR = mesh.nodes.coord(connR,:);

    phiR = phi(connR);

end


%% --------------------------------------------------------------
% Edge geometry
% --------------------------------------------------------------

localEdgeNodes = [ ...
    1 2;
    2 3;
    3 4;
    4 1];

n1 = localEdgeNodes(locL,1);
n2 = localEdgeNodes(locL,2);

x1 = ptsL(n1,:);
x2 = ptsL(n2,:);

dx = x2-x1;

edgeLength = norm(dx);

if edgeLength <= 0
    error('Zero-length edge %d.',edgeID);
end

edgeJac = edgeLength/2;


%% --------------------------------------------------------------
% Normal
% --------------------------------------------------------------

t = dx/edgeLength;

n = [t(2);-t(1)];


%% --------------------------------------------------------------
% Quadrature
% --------------------------------------------------------------

W = quad.edgeW;
Q = quad.edgeQ;

nGP = quad.nEdgeGP;


%% --------------------------------------------------------------
% Accumulation
% --------------------------------------------------------------

jumpIntegral = 0;


%% ==============================================================
% Edge integration
% ==============================================================

for q = 1:nGP

    s = Q(q);

    w = W(q);


    %% ----------------------------------------------------------
    % LEFT
    % ----------------------------------------------------------

    [xiL,etaL] = ...
        refEdgeCoords(locL,s);

    [N_L,dNdxiL] = ...
        lagrange_basis( ...
        mesh.meta.elemType, ...
        [xiL,etaL]);

    J_L = ptsL'*dNdxiL;

    dNdxL = dNdxiL/J_L;

    gradPhiL = dNdxL'*phiL;

    fluxL = ...
        -material.Gc*material.l0 ...
        *(gradPhiL'*n);


    %% ----------------------------------------------------------
    % RIGHT
    % ----------------------------------------------------------

    if hasRight

        [xiR,etaR] = ...
            refEdgeCoords(locR,s);

        [N_R,dNdxiR] = ...
            lagrange_basis( ...
            mesh.meta.elemType, ...
            [xiR,etaR]);

        J_R = ptsR'*dNdxiR;

        dNdxR = dNdxiR/J_R;

        gradPhiR = dNdxR'*phiR;

        fluxR = ...
            -material.Gc*material.l0 ...
            *(gradPhiR'*n);

    else

        % Homogeneous Neumann boundary
        fluxR = 0;

    end


    %% ----------------------------------------------------------
    % Flux jump
    % ----------------------------------------------------------

    jump = fluxL-fluxR;


    %% ----------------------------------------------------------
    % Integral
    % ----------------------------------------------------------

    jumpIntegral = ...
        jumpIntegral ...
        + jump^2*w*edgeJac;

end


%% --------------------------------------------------------------
% Verfürth scaling
%
% h_F = |F| in 2D
% --------------------------------------------------------------

contribution = ...
    edgeLength*jumpIntegral;

end