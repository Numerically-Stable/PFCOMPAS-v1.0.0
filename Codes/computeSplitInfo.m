function split = computeSplitInfo(strainVoigt,material,model)
% ============================================================
% computeSplitInfo - Constitutive Tension-Compression Splits
%
% PHASE-FIELD PHYSICS RELEVANCE:
% Standard isotropic phase-field models degrade the elastic strain energy 
% uniformly. Physically, this is catastrophic for general stress states: 
% it allows cracks to propagate under pure hydrostatic compression and 
% permits cracked element faces to physically interpenetrate.
%
% To enforce tension-compression asymmetry, the intact strain energy 
% psi_0(eps) must be decomposed into active (tensile, psi+) and passive 
% (compressive, psi-) components. The degradation function g(phi) is then 
% applied strictly to psi+, while psi- remains intact.
%
% This function robustly evaluates this split for both the stress tensor 
% and the 4th-order tangent stiffness matrix, supporting the computationally 
% efficient Amor (volumetric-deviatoric) split and the highly accurate 
% Miehe Spectral (principal strain) split.
% ============================================================
E  = material.E;
nu = material.nu;
dim = model.dim;

%% -----------------------------------------------------------
% Lamé constants
% ------------------------------------------------------------
mu     = E/(2*(1+nu));
lambda = E*nu/((1+nu)*(1-2*nu));

%% -----------------------------------------------------------
% Build full isotropic elastic stiffness matrix (C0)
% ------------------------------------------------------------
C = matStiff(E,nu,dim);

%% -----------------------------------------------------------
% Convert engineering strain vector to 2nd-order tensor
% ------------------------------------------------------------
% Note the division of the shear components by 2. This maps the 
% engineering shear strains (gamma) from the Voigt vector back to 
% the rigorous tensorial shear strains (epsilon_ij).
if dim == 2
    eps = [ ...
        strainVoigt(1),       0.5*strainVoigt(3), 0;
        0.5*strainVoigt(3),  strainVoigt(2),     0;
        0,                  0,                   0 ];
elseif dim == 3
    eps = [ ...
        strainVoigt(1),       0.5*strainVoigt(4), 0.5*strainVoigt(6);
        0.5*strainVoigt(4),  strainVoigt(2),     0.5*strainVoigt(5);
        0.5*strainVoigt(6),  0.5*strainVoigt(5), strainVoigt(3) ];
else
    error('Unsupported spatial dimension.');
end

%% -----------------------------------------------------------
% Full (Intact) elastic response
% ------------------------------------------------------------
trEps = trace(eps);
sigma = lambda*trEps*eye(3) + 2*mu*eps;
psi   = 0.5*lambda*trEps^2 + mu*sum(eps(:).^2);

%% -----------------------------------------------------------
% Pre-allocate outputs
% ------------------------------------------------------------
epsPlus  = zeros(3,3); epsMinus = zeros(3,3);
sigmaPlus  = zeros(3,3); sigmaMinus = zeros(3,3);
psiPlus  = 0; psiMinus = 0;
Cplus  = zeros(size(C)); Cminus = zeros(size(C));

%% ===========================================================
% 1. NO SPLIT (Isotropic)
% ============================================================
if strcmpi(model.split,'none')
    epsPlus   = eps;
    sigmaPlus = sigma;
    psiPlus   = psi;
    Cplus = C;

%% ===========================================================
% 2. AMOR SPLIT (Volumetric-Deviatoric)
% ============================================================
% Divides the energy into volumetric expansion (which can cause fracture) 
% and volumetric compression (which cannot). The deviatoric shear strain 
% energy is assumed to always contribute to fracture.
elseif strcmpi(model.split,'amor')
    I = eye(3);
    epsVol = (trEps/3)*I;
    epsDev = eps - epsVol;
    Kbulk = lambda + 2*mu/3;
    
    % --------------------------------------------------------
    % Algorithmic Regularization (Smoothing)
    % --------------------------------------------------------
    % The standard Macaulay bracket <x>+ = (x + |x|)/2 is not 
    % differentiable at x=0. To maintain Newton-Raphson stability, 
    % we replace |x| with sqrt(x^2 + eta^2), yielding a C1-continuous 
    % smooth approximation.
    eta = getSplitEta(eps,model);
    trPlus = 0.5*(trEps + sqrt(trEps^2 + eta^2));
    trMinus = trEps - trPlus;
    
    % Derivatives of the smooth Macaulay bracket w.r.t the trace
    dtrPlus = 0.5*(1 + trEps / sqrt(trEps^2 + eta^2));
    dtrMinus = 1 - dtrPlus;
    
    % --------------------------------------------------------
    % Energy and Stress Split
    % --------------------------------------------------------
    psiPlus  = 0.5*Kbulk*trPlus^2 + mu*sum(epsDev(:).^2);
    psiMinus = 0.5*Kbulk*trMinus^2;
    
    epsPlus  = epsDev + (trPlus/3)*I;
    epsMinus = (trMinus/3)*I;
    
    sigmaPlus  = Kbulk*trPlus*I + 2*mu*epsDev;
    sigmaMinus = Kbulk*trMinus*I;
    
    % --------------------------------------------------------
    % Consistent Amor Tangent
    % --------------------------------------------------------
    Cplus = buildAmorTangent(Kbulk,mu,dtrPlus,dim);
    Cminus = buildVolumetricTangent(Kbulk,dtrMinus,dim);

%% ===========================================================
% 3. MIEHE SPECTRAL SPLIT
% ============================================================
% The most physically rigorous split. Decomposes the strain tensor into 
% principal directions. Only positive principal strains (pure tension) 
% are degraded. Strictly prevents crack face interpenetration.
elseif strcmpi(model.split,'spectral')
    % --------------------------------------------------------
    % Eigen-decomposition
    % --------------------------------------------------------
    [V,D] = eig(eps);
    principalStrain = diag(D);
    
    % --------------------------------------------------------
    % Smooth Macaulay Bracket for Principal Strains
    % --------------------------------------------------------
    eta = getSplitEta(eps,model);
    lambdaPlus = 0.5*(principalStrain + sqrt(principalStrain.^2 + eta^2));
    lambdaMinus = principalStrain - lambdaPlus;
    
    % --------------------------------------------------------
    % Tensor Reconstruction
    % --------------------------------------------------------
    for a = 1:3
        v = V(:,a);
        P = v*v'; % Eigen-projection tensor (n_a x n_a)
        epsPlus = epsPlus + lambdaPlus(a)*P;
        epsMinus = epsMinus + lambdaMinus(a)*P;
    end
    
    % Volumetric parts for the Lambda term in the energy functional
    trEps = trace(eps);
    trPlus = 0.5*(trEps + sqrt(trEps^2 + eta^2));
    trMinus = trEps - trPlus;
    
    % --------------------------------------------------------
    % Energy and Stress Split
    % --------------------------------------------------------
    psiPlus  = 0.5*lambda*trPlus^2 + mu*sum(epsPlus(:).^2);
    psiMinus = 0.5*lambda*trMinus^2 + mu*sum(epsMinus(:).^2);
    
    sigmaPlus  = lambda*trPlus*eye(3) + 2*mu*epsPlus;
    sigmaMinus = lambda*trMinus*eye(3) + 2*mu*epsMinus;
    
    % --------------------------------------------------------
    % Consistent Spectral Tangents
    % --------------------------------------------------------
    % The derivative d(sigma_plus)/d(eps) requires complex 4th-order tensor 
    % calculus, taking into account the rotation of the principal axes.
    Cplus = buildSpectralTangent(V,principalStrain,lambdaPlus,lambda,mu,eta,'plus',dim);
    Cminus = buildSpectralTangent(V,principalStrain,lambdaMinus,lambda,mu,eta,'minus',dim);

else
    error('Unknown split: %s',model.split);
end

%% -----------------------------------------------------------
% Assembly of returned Split Struct
% ------------------------------------------------------------
Csplit = Cplus + Cminus;

split.eps = eps; split.epsPlus = epsPlus; split.epsMinus = epsMinus;
split.sigma = sigma; split.sigmaPlus = sigmaPlus; split.sigmaMinus = sigmaMinus;
split.psi = psi; split.psiPlus = psiPlus; split.psiMinus = psiMinus;
split.C = C; split.Cplus = Cplus; split.Cminus = Cminus; split.Csplit = Csplit;

end