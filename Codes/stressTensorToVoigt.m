function sigmaVoigt = stressTensorToVoigt(sigma,dim)
% ============================================================
% stressTensorToVoigt - Map Cauchy Stress to Voigt Vector
%
% COMPUTATIONAL MECHANICS RELEVANCE (WORK-CONJUGACY):
% To maintain thermodynamic consistency, the internal virtual work evaluated 
% using Voigt vectors must identically match the double contraction of the 
% physical tensors:
%       W_int = sigma : epsilon = sigma_voigt^T * epsilon_voigt
%
% Tensorial expansion (2D):
%       sigma : epsilon = sigma_xx*eps_xx + sigma_yy*eps_yy + 2*sigma_xy*eps_xy
%
% Because the B-matrix outputs the ENGINEERING shear strain (gamma_xy = 2*eps_xy),
% the Voigt expansion evaluates as:
%       sigma_voigt^T * epsilon_voigt = sigma_xx*eps_xx + sigma_yy*eps_yy + tau_xy*gamma_xy
%
% Therefore, to make these two equations equal, the shear stress component 
% in the Voigt vector MUST be the raw Cauchy shear stress (tau_xy = sigma_xy). 
% There is NO factor of 2 for shear stresses.
% ============================================================

if dim == 2
    % --------------------------------------------------------
    % 2D Voigt Stress Vector
    % Ordering: [sigma_xx; sigma_yy; tau_xy]
    % --------------------------------------------------------
    sigmaVoigt = [
        sigma(1,1);
        sigma(2,2);
        sigma(1,2)
    ];
    
elseif dim == 3
    % --------------------------------------------------------
    % 3D Voigt Stress Vector
    % Ordering: [sigma_xx; sigma_yy; sigma_zz; tau_xy; tau_yz; tau_xz]
    % --------------------------------------------------------
    sigmaVoigt = [
        sigma(1,1);
        sigma(2,2);
        sigma(3,3);
        sigma(1,2);  % tau_xy
        sigma(2,3);  % tau_yz
        sigma(1,3)   % tau_xz
    ];
else
    error('Unsupported dimension dim = %d. Use dim = 2 or 3.', dim);
end

end