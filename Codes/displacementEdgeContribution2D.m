function [eL,eR,contribution] = ...
    displacementEdgeContribution2D( ...
    mesh,edgeID,u,phi,material,model,quad)
% ==============================================================
% DISPLACEMENT FLUX JUMP ON ONE 2D EDGE
%
% Returns one scalar contribution.
% ==============================================================


%% --------------------------------------------------------------
% Connectivity
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

uL = u(mesh.elem.dofU(eL,:));

phiL = phi(connL);


%% --------------------------------------------------------------
% Right element
% --------------------------------------------------------------

if hasRight

    locR = mesh.edge.conn(edgeID,4);

    connR = mesh.elem.conn(eR,:);

    ptsR = mesh.nodes.coord(connR,:);

    uR = u(mesh.elem.dofU(eR,:));

    phiR = phi(connR);

end


%% --------------------------------------------------------------
% Geometry
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
% Integration
% ==============================================================

for q = 1:nGP

    s = Q(q);

    w = W(q);


    %% ==========================================================
    % LEFT
    % ==========================================================

    [xiL,etaL] = ...
        refEdgeCoords(locL,s);

    [N_L,dNdxiL] = ...
        lagrange_basis( ...
        mesh.meta.elemType, ...
        [xiL,etaL]);

    J_L = ptsL'*dNdxiL;

    dNdxL = dNdxiL/J_L;


    % ----------------------------------------------------------
    % B matrix
    % ----------------------------------------------------------

    B_L = displacementBMatrix( ...
        dNdxL, ...
        model.dim, ...
        mesh.elem.nen);


    % ----------------------------------------------------------
    % Strain
    % ----------------------------------------------------------

    strainL = B_L*uL;


    % ----------------------------------------------------------
    % Phase field
    % ----------------------------------------------------------

    phiL_q = N_L'*phiL;


    % ----------------------------------------------------------
    % Constitutive response
    % ----------------------------------------------------------

    sigmaL = constitutiveStress( ...
        strainL, ...
        phiL_q, ...
        material, ...
        model);


    % ----------------------------------------------------------
    % Traction
    % ----------------------------------------------------------
    
    sigmaMat = [sigmaL(1) sigmaL(3);
                    sigmaL(3) sigmaL(2)];
    tractionL = sigmaMat*n;


    %% ==========================================================
    % RIGHT
    % ==========================================================

    if hasRight

        [xiR,etaR] = ...
            refEdgeCoords(locR,s);

        [N_R,dNdxiR] = ...
            lagrange_basis( ...
            mesh.meta.elemType, ...
            [xiR,etaR]);

        J_R = ptsR'*dNdxiR;

        dNdxR = dNdxiR/J_R;


        B_R = displacementBMatrix( ...
            dNdxR, ...
            model.dim, ...
            mesh.elem.nen);


        strainR = B_R*uR;

        phiR_q = N_R'*phiR;


        sigmaR = constitutiveStress( ...
            strainR, ...
            phiR_q, ...
            material, ...
            model);

        
        sigmaMat = [sigmaR(1) sigmaR(3);
                    sigmaR(3) sigmaR(2)];
        tractionR = sigmaMat*n;

    else

        tractionR = zeros(2,1);

    end


    %% ----------------------------------------------------------
    % Jump
    % ----------------------------------------------------------

    jump = tractionL-tractionR;


    %% ----------------------------------------------------------
    % Integration
    % ----------------------------------------------------------

    jumpIntegral = ...
        jumpIntegral ...
        + dot(jump,jump)*w*edgeJac;

end


%% --------------------------------------------------------------
% Verfürth scaling
% --------------------------------------------------------------

contribution = ...
    edgeLength*jumpIntegral;

end