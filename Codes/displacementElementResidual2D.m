function eta2 = displacementElementResidual2D( ...
    mesh,e,u,phi,material,model,quad)
% ==============================================================
% Element residual contribution.
%
% Current formulation does not reconstruct div(sigma).
% Therefore:
%
%       eta_K^2 = 0
%
% The jump residual remains the active contribution.
% ==============================================================

eta2 = 0.0;

end