function H = voigtBasisTensor(j,dim)
% ============================================================
% voigtBasisTensor - Generate 2nd-Order Symmetric Unit Basis Tensors
%
% COMPUTATIONAL MECHANICS RELEVANCE:
% This function generates the symmetric 2nd-order strain tensor H 
% corresponding to a unit perturbation in the j-th Voigt direction.
% It is utilized by the tangent extraction routines (like buildAmorTangent) 
% to programmatically build the 4th-order constitutive tensor.
%
% ENGINEERING SHEAR STRAIN COUPLING:
% A unit perturbation in normal strain (e.g., epsilon_xx = 1) translates 
% directly to H(1,1) = 1.
% 
% However, the Voigt vector uses ENGINEERING shear strains (gamma = 2*eps).
% Therefore, a unit perturbation in the Voigt engineering shear strain 
% (e.g., gamma_xy = 1) means the actual physical tensorial shear strain 
% must be perturbed by exactly half that amount: 
%       epsilon_xy = epsilon_yx = 0.5
%
% This function explicitly handles this 0.5 mapping, guaranteeing that 
% the resulting tangent moduli perfectly match the energy expectations 
% of the B-matrix.
% ============================================================

H = zeros(3,3);

if dim == 2
    switch j
        case 1
            % Unit perturbation in normal x-strain (epsilon_xx)
            H(1,1) = 1;
        case 2
            % Unit perturbation in normal y-strain (epsilon_yy)
            H(2,2) = 1;
        case 3
            % Unit perturbation in engineering shear strain (gamma_xy)
            % epsilon_xy = epsilon_yx = 0.5 * gamma_xy
            H(1,2) = 0.5;
            H(2,1) = 0.5;
        otherwise
            error('Invalid 2D Voigt index.');
    end
    
elseif dim == 3
    switch j
        case 1
            H(1,1) = 1; % epsilon_xx
        case 2
            H(2,2) = 1; % epsilon_yy
        case 3
            H(3,3) = 1; % epsilon_zz
        case 4
            % gamma_xy
            H(1,2) = 0.5;
            H(2,1) = 0.5;
        case 5
            % gamma_yz
            H(2,3) = 0.5;
            H(3,2) = 0.5;
        case 6
            % gamma_xz
            H(1,3) = 0.5;
            H(3,1) = 0.5;
        otherwise
            error('Invalid 3D Voigt index.');
    end
else
    error('Unsupported dimension.');
end
end