function [D] = matStiff(E, nu, dim)
% ============================================================
% matStiff - Pristine Isotropic Material Stiffness Matrix (Voigt)
%
% COMPUTATIONAL MECHANICS RELEVANCE:
% This function formulates the baseline, undamaged 4th-order isotropic 
% elasticity tensor (C_0) mapped into a 2D matrix (D) via Voigt notation.
% 
% For phase-field formulations utilizing the 'none' (isotropic) split, 
% this exact matrix is simply multiplied by the degradation function g(phi). 
% For asymmetric splits (Amor/Spectral), this matrix provides the baseline 
% elastic bounds.
%
% IMPORTANT - ENERGY EQUIVALENCE:
% As established in the displacementBMatrix function, the framework utilizes 
% engineering shear strains (gamma_xy = 2*epsilon_xy). To preserve 
% exact thermodynamic energy equivalence such that:
%     1/2 * epsilon : C_0 : epsilon == 1/2 * epsilon_voigt^T * D * epsilon_voigt
% the diagonal shear components of D must strictly equal the shear 
% modulus (G or mu), NOT 2*G or G/2.
% ============================================================

mu0 = E / (2 * (1 + nu));           % Shear modulus (G or mu)
K0 = E / (3 * (1 - 2 * nu));        % Bulk modulus (K)
lambda0 = K0 - 2 * mu0 / 3;         % Lamé's first parameter (lambda)

if dim == 3  
    % --------------------------------------------------------
    % Full 3D Elasticity Matrix
    % Voigt Ordering: [xx, yy, zz, xy, yz, xz]
    % --------------------------------------------------------
    D = zeros(6, 6);
    for i = 1:3
        for j = 1:3
            D(j, i) = lambda0;
        end
        % Diagonal normal terms: lambda + 2*mu
        D(i, i) = 2 * mu0 + lambda0;
    end
    for i = 4:6
        % Diagonal shear terms: mu
        D(i, i) = mu0;
    end
    
elseif dim == 2  
    % --------------------------------------------------------
    % 2D Plane Strain Matrix
    % Voigt Ordering: [xx, yy, xy]
    %
    % NOTE: Phase-field fracture in 2D is almost universally modeled 
    % under plane strain conditions rather than plane stress. Plane stress 
    % implies out-of-plane Poisson contraction (necking), which physically 
    % interferes with the interpretation of the characteristic length scale l0 
    % and the critical fracture energy Gc.
    % --------------------------------------------------------
    D = zeros(3, 3);
    factor = E / ((1 + nu) * (1 - 2 * nu));
    
    % Normal stress/strain coupling
    D(1,1) = factor * (1 - nu);
    D(1,2) = factor * nu;
    D(2,1) = D(1,2);
    D(2,2) = factor * (1 - nu);
    
    % Shear coupling
    G = E / (2 * (1 + nu));   % Shear modulus
    D(3,3) = G;               % Strictly G (not 2G) for engineering shear
else
    error('Invalid dim. Must be 2 or 3.');
end

end