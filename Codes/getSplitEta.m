function eta = getSplitEta(eps,model)
% ============================================================
% getSplitEta - Algorithmic Smoothing Parameter for Energy Splits
%
% PHASE-FIELD HPC RELEVANCE:
% To enforce tension-compression asymmetry (via Amor or Spectral splits), 
% the solver must isolate the positive (tensile) and negative (compressive) 
% components of the strain tensor or its trace. 
%
% Mathematically, this is done using the Macaulay bracket: 
%       <x>_pm = 0.5 * (x +/- |x|)
%
% The absolute value function |x|, however, is non-differentiable at x=0. 
% This sharp "kink" implies that the strain energy density is only C^0 
% continuous. Consequently, the stress tensor (first derivative) is 
% discontinuous, and the tangent stiffness matrix (second derivative) 
% is mathematically undefined at the transition point between tension 
% and compression. If a Gauss point lands exactly on (or oscillates across) 
% this boundary, the Newton-Raphson solver will catastrophically fail.
%
% To guarantee asymptotic quadratic convergence, we replace |x| with a 
% smooth, C1-continuous hyperbola: sqrt(x^2 + eta^2).
% This function calculates that eta parameter.
% ============================================================

%% -----------------------------------------------------------
% Default Regularization Base (eta0)
% ------------------------------------------------------------
% A vanishingly small base parameter (10^-12) ensures the physical 
% accuracy of the split is not compromised by over-smoothing.
if isfield(model,'splitEta')
    eta0 = model.splitEta;
else
    eta0 = 1e-12;
end

%% -----------------------------------------------------------
% Dynamic Strain Scaling
% ------------------------------------------------------------
% If the local strain magnitude becomes exceptionally large, a fixed 
% eta = 1e-12 might fall below the relative floating-point precision 
% limits of the machine (eps ~ 2.22e-16), effectively un-smoothing the 
% bracket and causing NaN errors in the tangent.
%
% Conversely, we clamp the minimum scale to 1.0 to prevent eta from 
% dropping to zero in entirely unstrained regions.
strainScale = max(abs(eps(:)));
strainScale = max(1,strainScale);

%% -----------------------------------------------------------
% Final Regularization Parameter
% ------------------------------------------------------------
% The regularization is dynamically scaled with the local deformation state.
eta = eta0 * strainScale;

end