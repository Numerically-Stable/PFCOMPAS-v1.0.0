function [Eelastic,Efracture] = ...
    elementAT2Energy( ...
    mesh,e,u,phi,material,model,quad)
% ============================================================
% Element AT2 energy contribution.
%
% Called only during postprocessing.
% ============================================================

dim = mesh.meta.dim;

conn   = mesh.elem.conn(e,:);
coords = mesh.nodes.coord(conn,:);

dofU = mesh.elem.dofU(e,:);

u_e   = u(dofU);
phi_e = phi(conn);

nen = mesh.elem.nen;

W = quad.W;
Q = quad.Q;

ngp = numel(W);

Gc  = material.Gc;
ell = material.l0;


%% ============================================================
% Initialize
% ============================================================

Eelastic  = 0;
Efracture = 0;


%% ============================================================
% Gauss loop
% ============================================================

for q = 1:ngp

    %% --------------------------------------------------------
    % Shape functions
    % ---------------------------------------------------------

    [N,dNdxi] = ...
        lagrange_basis( ...
        mesh.meta.elemType,Q(q,:));

    %% --------------------------------------------------------
    % Geometry
    % ---------------------------------------------------------

    J = coords' * dNdxi;

    detJ = det(J);

    if detJ <= 0

        error( ...
            'Non-positive Jacobian: element %d, GP %d.', ...
            e,q);

    end

    dNdx = dNdxi / J;

    %% --------------------------------------------------------
    % Displacement strain
    % ---------------------------------------------------------

    B = displacementBMatrix( ...
        dNdx,dim,nen);

    strain = B*u_e;

    %% --------------------------------------------------------
    % Phase field
    % ---------------------------------------------------------

    phi_gp = N'*phi_e;

    gradPhi = dNdx'*phi_e;

    %% --------------------------------------------------------
    % Degradation
    % ---------------------------------------------------------

    [g,~,~] = ...
        degradationFunction( ...
        phi_gp,material);

    %% --------------------------------------------------------
    % Constitutive split
    % ---------------------------------------------------------

    split = computeSplitInfo( ...
        strain,material,model);

    %% --------------------------------------------------------
    % Elastic energy density
    % ---------------------------------------------------------

    psi = ...
        g*split.psiPlus ...
        + split.psiMinus;

    %% --------------------------------------------------------
    % Fracture energy density
    % ---------------------------------------------------------

    gamma = ...
        1/(2*ell)*phi_gp^2 ...
        + ell/2*(gradPhi'*gradPhi);

    fractureDensity = ...
        Gc*gamma;

    %% --------------------------------------------------------
    % Integration
    % ---------------------------------------------------------

    weight = W(q)*detJ;

    Eelastic = ...
        Eelastic + psi*weight;

    Efracture = ...
        Efracture + fractureDensity*weight;

end

end