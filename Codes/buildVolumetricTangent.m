function Cvoigt = buildVolumetricTangent(K,dtr,dim)
% ============================================================
% buildVolumetricTangent - Tangent Moduli for Volumetric Compression
%
% PHASE-FIELD HPC RELEVANCE:
% In the Amor volumetric-deviatoric split, the compressive (passive) 
% strain energy is assumed to be purely volumetric. Consequently, the 
% passive stress tensor reduces to a purely dilatational state:
%       sigma_minus = K * <tr(eps)>_minus * I
%
% To maintain the asymptotic quadratic convergence of the Newton-Raphson 
% solver across compressive zones, the exact analytical derivative of this 
% passive stress must be evaluated to form the C_minus tangent matrix.
%
% Using the same highly robust kinematic perturbation technique as 
% buildAmorTangent, this function sequentially feeds unit Voigt strain 
% tensors (H) into the linearized stress equation to extract the exact 
% tangent column by column, bypassing explicit 4th-order algebraic unrolling.
% ============================================================

%% Number of Voigt components
if dim == 2
    nVoigt = 3;
elseif dim == 3
    nVoigt = 6;
else
    error('Unsupported dimension dim = %d.', dim);
end

Cvoigt = zeros(nVoigt,nVoigt);
I = eye(3);

%% ============================================================
% Kinematic Perturbation Loop
% ============================================================
for j = 1:nVoigt
    % --------------------------------------------------------
    % Engineering strain perturbation (Basis Tensor)
    % --------------------------------------------------------
    H = voigtBasisTensor(j,dim);
    trH = trace(H);
    
    % --------------------------------------------------------
    % Volumetric Stress Perturbation (Analytical Derivative)
    % --------------------------------------------------------
    % Evaluates the exact change in the Amor passive stress tensor 
    % resulting from the unit basis strain H.
    % 
    % Delta_sigma_minus = K * [ d(tr_minus)/d(tr_eps) * tr(H) ] * I 
    dsigma = K*dtr*trH*I;
    
    % --------------------------------------------------------
    % Convert tensor stress to Voigt notation (Column Mapping)
    % --------------------------------------------------------
    % The symmetric stress perturbation is flattened into a Voigt vector 
    % and assigned to the j-th column of the compressive tangent matrix.
    Cvoigt(:,j) = stressTensorToVoigt(dsigma,dim);
end
end