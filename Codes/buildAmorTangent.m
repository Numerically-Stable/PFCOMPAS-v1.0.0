function Cvoigt = buildAmorTangent(K,mu,dtr,dim)
% ============================================================
% buildAmorTangent - Tangent Moduli for Volumetric-Deviatoric Split
%
% PHASE-FIELD HPC RELEVANCE:
% To achieve quadratic convergence in the displacement Newton-Raphson solver, 
% we require the exact algorithmic tangent stiffness: C = d(sigma)/d(eps).
% 
% For the Amor split, the active (plus) stress is:
%       sigma_plus = K * <tr(eps)>_plus * I + 2 * mu * eps_dev
%
% Computing the 4th-order derivative analytically and mapping it to Voigt 
% notation is highly error-prone. Instead, this function uses a programmatic 
% kinematic perturbation approach. By sequentially feeding unit "basis strains" 
% (H) into the analytical derivative equation, the 2D Voigt tangent matrix 
% is robustly populated column by column.
%
% INPUT
%   K   : Bulk modulus
%   mu  : Shear modulus
%   dtr : Derivative of the smoothed positive/negative volumetric strain trace
%   dim : Spatial dimension (2 or 3)
%
% OUTPUT
%   Cvoigt : Consistent tangent matrix mapped to engineering Voigt notation
%            (3x3 for 2D, 6x6 for 3D)
% ============================================================

%% Number of Voigt components
if dim == 2
    nVoigt = 3;
elseif dim == 3
    nVoigt = 6;
else
    error('Unsupported dimension dim = %d.', dim);
end

%% Initialize Tangent Matrix and Identity Tensor
Cvoigt = zeros(nVoigt,nVoigt);
I = eye(3);

%% ============================================================
% Kinematic Perturbation Loop
% ============================================================
% Iterate through each discrete Voigt direction.
for j = 1:nVoigt
    
    % --------------------------------------------------------
    % Engineering strain perturbation (Basis Tensor)
    % --------------------------------------------------------
    % Extracts a 3x3 symmetric tensor H corresponding to a unit strain 
    % in the j-th Voigt direction (e.g., if j=1, H has a 1 at (1,1); 
    % if j=3 in 2D, H represents a pure shear state).
    H = voigtBasisTensor(j,dim);
    trH = trace(H);
    
    % --------------------------------------------------------
    % Deviatoric perturbation
    % --------------------------------------------------------
    Hdev = H - (trH/3)*I;
    
    % --------------------------------------------------------
    % Stress perturbation (Analytical Derivative)
    % --------------------------------------------------------
    % Evaluates the exact change in the Amor active stress tensor 
    % resulting from the unit basis strain H.
    % 
    % Delta_sigma_plus = K * [ d(tr_plus)/d(tr_eps) * tr(H) ] * I 
    %                    + 2 * mu * H_dev
    dsigma = K*dtr*trH*I + 2*mu*Hdev;
    
    % --------------------------------------------------------
    % Convert tensor stress to Voigt notation (Column Mapping)
    % --------------------------------------------------------
    % The resulting symmetric stress perturbation tensor is flattened 
    % into a Voigt vector and strictly assigned to the j-th column 
    % of the global tangent moduli matrix.
    Cvoigt(:,j) = stressTensorToVoigt(dsigma,dim);
end

end