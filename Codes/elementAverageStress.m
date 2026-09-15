function sigma_avg = elementAverageStress( ...
    mesh,e,u,phi, ...
    material,model,quad)

% ==============================================================
% ELEMENT AVERAGE STRESS
%
% Computes the volume/area averaged Cauchy stress for element e:
%
%       sigma_avg = (1/|Omega_e|) * integral_Omega sigma dOmega
%
% 2D output:
%
%       sigma_avg = [sigma_xx;
%                    sigma_yy;
%                    sigma_xy]
%
% 3D output:
%
%       sigma_avg = [sigma_xx;
%                    sigma_yy;
%                    sigma_zz;
%                    sigma_xy;
%                    sigma_yz;
%                    sigma_xz]
%
% The constitutive stress is evaluated using:
%
%       1. displacement u
%       2. phase field phi
%       3. material
%       4. model.split
%
% ===============================================================


%% ==============================================================
% BASIC INFORMATION
% ==============================================================

dim = model.dim;
nen = mesh.elem.nen;

ngp = numel(quad.W);


%% ==============================================================
% ELEMENT DATA
% ==============================================================

conn = mesh.elem.conn(e,:);

pts = mesh.nodes.coord(conn,:);

dofU = mesh.elem.dofU(e,:);

u_e = u(dofU);

phi_e = phi(conn);


%% ==============================================================
% INITIALIZE
% ==============================================================

if dim == 2

    nStress = 3;

elseif dim == 3

    nStress = 6;

else

    error( ...
        'Unsupported spatial dimension: %d', ...
        dim);

end


sigma_integral = zeros(nStress,1);

volume = 0;


%% ==============================================================
% GAUSS INTEGRATION
% ==============================================================

for q = 1:ngp

    pt = quad.Q(q,:);
    wt = quad.W(q);


    %% ----------------------------------------------------------
    % Shape functions
    %% ----------------------------------------------------------

    [N,dNdxi] = ...
        lagrange_basis( ...
        mesh.meta.elemType, ...
        pt);


    %% ----------------------------------------------------------
    % Jacobian
    %% ----------------------------------------------------------

    J = pts' * dNdxi;

    detJ = det(J);


    if detJ <= 0

        error( ...
            ['Non-positive Jacobian in element %d ' ...
             'at Gauss point %d: detJ = %.6e'], ...
            e,q,detJ);

    end


    %% ----------------------------------------------------------
    % Spatial derivatives
    %% ----------------------------------------------------------

    dNdx = dNdxi / J;


    %% ----------------------------------------------------------
    % Strain-displacement matrix
    %% ----------------------------------------------------------

    B = displacementBMatrix( ...
        dNdx, ...
        dim, ...
        nen);


    %% ----------------------------------------------------------
    % Strain
    %% ----------------------------------------------------------

    strain = B * u_e;


    %% ----------------------------------------------------------
    % Phase field
    %% ----------------------------------------------------------

    phi_gp = N' * phi_e;


    %% ----------------------------------------------------------
    % Degradation
    %% ----------------------------------------------------------

    [g,~,~] = ...
        degradationFunction( ...
        phi_gp, ...
        model, ...
        material);


    %% ----------------------------------------------------------
    % Constitutive split
    %% ----------------------------------------------------------

    split = computeSplitInfo( ...
        strain, ...
        material, ...
        model);


    %% ----------------------------------------------------------
    % Stress
    %% ----------------------------------------------------------

    switch lower(model.split)

        case 'none'

            % Undegraded constitutive stress multiplied by g.
            sigma = g * split.sigma;


        case {'amor','spectral'}

            % Tensile part degraded, compressive part retained.
            sigma = ...
                g * split.sigmaPlus ...
                + split.sigmaMinus;


        otherwise

            error( ...
                'Unknown stress split: %s', ...
                model.split);

    end


    %% ----------------------------------------------------------
    % Convert tensor to FEM Voigt ordering
    %% ----------------------------------------------------------

    sigmaVoigt = ...
        stressTensorToVoigt( ...
        sigma, ...
        dim);


    %% ----------------------------------------------------------
    % Accumulate integral
    %% ----------------------------------------------------------

    sigma_integral = ...
        sigma_integral ...
        + sigmaVoigt * wt * detJ;


    volume = ...
        volume + wt * detJ;

end


%% ==============================================================
% VOLUME AVERAGE
% ==============================================================

sigma_avg = ...
    sigma_integral / volume;

end