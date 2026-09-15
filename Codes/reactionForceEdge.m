function R_e = reactionForceEdge( ...
    mesh,edgeID,u,phi,material,model, ...
    W,Q,nGP,boundaryName,dir)
% ==============================================================
% SINGLE BOUNDARY EDGE REACTION CONTRIBUTION
%
% Returns only one scalar:
%
%       R_e = integral_edge (sigma*n)_dir dGamma
%
% Designed as an independent kernel so that the outer edge loop
% can later be replaced by PARFOR.
% ==============================================================


%% --------------------------------------------------------------
% Element / edge information
% --------------------------------------------------------------

e   = mesh.edge.conn(edgeID,1);
loc = mesh.edge.conn(edgeID,2);

conn = mesh.elem.conn(e,:);

coords = mesh.nodes.coord(conn,:);

dofU = mesh.elem.dofU(e,:);

u_e   = u(dofU);
phi_e = phi(conn);

nen = mesh.elem.nen;


%% --------------------------------------------------------------
% Local edge node numbering
% --------------------------------------------------------------

localEdgeNodes = [ ...
    1 2;
    2 3;
    3 4;
    4 1];


edgeNodes = conn(localEdgeNodes(loc,:));


%% --------------------------------------------------------------
% Edge geometry
% --------------------------------------------------------------

x1 = mesh.nodes.coord(edgeNodes(1),:);
x2 = mesh.nodes.coord(edgeNodes(2),:);

dx = x2 - x1;

edgeLength = norm(dx);

if edgeLength <= 0

    error( ...
        'Zero-length boundary edge %d.',edgeID);

end

edgeJac = 0.5*edgeLength;


%% --------------------------------------------------------------
% Tangent
% --------------------------------------------------------------

t = dx/edgeLength;


%% --------------------------------------------------------------
% Outward normal
%
% Initial normal for counter-clockwise edge orientation.
% Boundary name is used to enforce the expected outward
% direction.
% --------------------------------------------------------------

n = [t(2);-t(1)];


switch lower(boundaryName)

    case 'top'

        if n(2) < 0
            n = -n;
        end

    case 'bottom'

        if n(2) > 0
            n = -n;
        end

    case 'right'

        if n(1) < 0
            n = -n;
        end

    case 'left'

        if n(1) > 0
            n = -n;
        end

    otherwise

        % Keep the geometrically constructed normal.
        %
        % For arbitrary boundary names, the mesh orientation
        % must therefore be consistent.

end


%% --------------------------------------------------------------
% Local scalar accumulation
% --------------------------------------------------------------

R_e = 0.0;


%% ==============================================================
% Edge Gauss integration
% ==============================================================

for q = 1:nGP

    %% ----------------------------------------------------------
    % Reference edge coordinate
    % ----------------------------------------------------------

    s = Q(q);


    %% ----------------------------------------------------------
    % Map 1D edge coordinate to element reference coordinates
    % ----------------------------------------------------------

    [xi,eta] = refEdgeCoords(loc,s);


    %% ----------------------------------------------------------
    % Shape functions
    % ----------------------------------------------------------

    [N,dNdxi] = ...
        lagrange_basis( ...
        mesh.meta.elemType, ...
        [xi,eta]);


    %% ----------------------------------------------------------
    % Element geometry
    % ----------------------------------------------------------

    J = coords' * dNdxi;

    detJ = det(J);

    if detJ <= 0

        error( ...
            ['Non-positive Jacobian: element %d, ' ...
             'edge %d, GP %d.'], ...
            e,edgeID,q);

    end


    %% ----------------------------------------------------------
    % Spatial derivatives
    % ----------------------------------------------------------

    dNdx = dNdxi/J;


    %% ----------------------------------------------------------
    % B matrix
    % ----------------------------------------------------------

    B = displacementBMatrix( ...
        dNdx, ...
        model.dim, ...
        nen);


    %% ----------------------------------------------------------
    % Strain
    % ----------------------------------------------------------

    strain = B*u_e;


    %% ----------------------------------------------------------
    % Phase field at GP
    % ----------------------------------------------------------

    phi_gp = N'*phi_e;


    %% ----------------------------------------------------------
    % Degradation
    % ----------------------------------------------------------

    [g,~,~] = ...
        degradationFunction( ...
        phi_gp, ...
        model, ...
        material);


    %% ----------------------------------------------------------
    % Constitutive split
    % ----------------------------------------------------------

    split = computeSplitInfo( ...
        strain, ...
        material, ...
        model);


    %% ----------------------------------------------------------
    % Stress
    % ----------------------------------------------------------

    switch lower(model.split)

        case 'none'

            sigma = ...
                g*split.sigma;


        case {'spectral','amor'}

            sigma = ...
                g*split.sigmaPlus ...
                + split.sigmaMinus;


        otherwise

            error( ...
                'Unknown split: %s', ...
                model.split);

    end


    %% ----------------------------------------------------------
    % Traction
    % ----------------------------------------------------------
    
    sigmaMat = [sigma(1) sigma(3);
                    sigma(3) sigma(2)];
    traction = sigmaMat*n;


    %% ----------------------------------------------------------
    % Boundary integration
    % ----------------------------------------------------------

    weight = W(q)*edgeJac;

    R_e = ...
        R_e ...
        + traction(dir)*weight;

end

end