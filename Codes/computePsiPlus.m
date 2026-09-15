function psiPlus = computePsiPlus(strainVoigt, material, model)
% ============================================================
% computePsiPlus
%
% Computes tensile elastic energy density psi_plus.
%
% Supported splits:
%   'none'
%   'spectral'
%   'amor'
%
% INPUT
%   strainVoigt : engineering strain vector
%                 2D -> [exx; eyy; gxy]
%                 3D -> [exx; eyy; ezz; gxy; gyz; gxz]
%
% OUTPUT
%   psiPlus     : tensile elastic energy density
%
% ============================================================

E  = material.E;
nu = material.nu;
dim = model.dim;

%% Convert engineering strain to 3x3 tensor

if dim == 2

    eps = [ ...
        strainVoigt(1),  0.5*strainVoigt(3), 0;
        0.5*strainVoigt(3), strainVoigt(2), 0;
        0,                  0,                0 ];

elseif dim == 3

    eps = [ ...
        strainVoigt(1),       0.5*strainVoigt(4), 0.5*strainVoigt(6);
        0.5*strainVoigt(4),  strainVoigt(2),      0.5*strainVoigt(5);
        0.5*strainVoigt(6),  0.5*strainVoigt(5),  strainVoigt(3) ];

else

    error('Unsupported spatial dimension.');

end

%% Lamé constants

mu     = E/(2*(1+nu));
lambda = E*nu/((1+nu)*(1-2*nu));

%% ==========================================================
% NO SPLIT
% ===========================================================

if strcmpi(model.split,'none')

    trEps = trace(eps);

    psiPlus = ...
        0.5*lambda*trEps^2 + ...
        mu*sum(eps(:).^2);

    return;

end

%% ==========================================================
% AMOR SPLIT
% ===========================================================

if strcmpi(model.split,'amor')

    trEps = trace(eps);

    epsDev = eps - (trEps/3)*eye(3);

    Kbulk = lambda + 2*mu/3;

    trPlus = max(trEps,0);

    psiPlus = ...
        0.5*Kbulk*trPlus^2 + ...
        mu*sum(epsDev(:).^2);

    return;

end

%% ==========================================================
% SPECTRAL SPLIT
% ===========================================================

if strcmpi(model.split,'spectral')

    [V,D] = eig(eps);

    principalStrain = diag(D);

    eta = 1e-12 * ...
        max(1,max(abs(principalStrain)));

    lambdaPlus = ...
        0.5*(principalStrain + ...
        sqrt(principalStrain.^2 + eta^2));

    epsPlus = zeros(3,3);

    for i = 1:3

        v = V(:,i);

        epsPlus = epsPlus + ...
            lambdaPlus(i)*(v*v');

    end

    % Correct volumetric tensile part
    trEps  = trace(eps);
    trPlus = max(trEps,0);

    psiPlus = ...
        0.5*lambda*trPlus^2 + ...
        mu*sum(epsPlus(:).^2);

    return;

end

%% Unknown split

error('Unknown phase-field split: %s',model.split);

end