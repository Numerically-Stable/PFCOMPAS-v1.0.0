function sigma = constitutiveStress( ...
    strain,phi_gp,material,model)
% ==============================================================
% CONSTITUTIVE STRESS
%
% Same constitutive path used by the displacement solver.
%
% Returns:
%
%   sigma : dim x dim stress tensor
% ==============================================================


%% --------------------------------------------------------------
% Degradation
% --------------------------------------------------------------

[g,~,~] = ...
    degradationFunction( ...
    phi_gp, ...
    model, ...
    material);


%% --------------------------------------------------------------
% Constitutive split
% --------------------------------------------------------------

split = computeSplitInfo( ...
    strain, ...
    material, ...
    model);


%% --------------------------------------------------------------
% Stress
% --------------------------------------------------------------

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
            'Unknown split: %s',model.split);

end

end