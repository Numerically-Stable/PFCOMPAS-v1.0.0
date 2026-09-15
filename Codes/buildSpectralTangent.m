function Cvoigt = buildSpectralTangent( ...
    V,lambdaEig,fEig,lambdaLam,mu,eta,branch,dim)
% ============================================================
% buildSpectralTangent - Consistent Tangent for Spectral Split
%
% PHASE-FIELD PHYSICS RELEVANCE (THE FRÉCHET DERIVATIVE):
% In the Miehe spectral split, the active/passive strain tensors are 
% defined by projecting the smoothed principal strains back onto the 
% principal directions (eigenvectors). 
%
% Taking the derivative of this tensor w.r.t the full strain tensor 
% is mathematically demanding. When a perturbation is applied to the 
% strain state, not only do the principal strain magnitudes (eigenvalues) 
% change, but the principal axes (eigenvectors) ROTATE. 
%
% To construct the exact consistent tangent, we must evaluate the 
% Fréchet derivative of this isotropic tensor function. This introduces 
% the concept of "divided differences" to account for the cross-coupling 
% between different principal directions during rotation.
%
% Using the same kinematic perturbation technique as Amor, we feed 
% unit basis strains (H) into this Fréchet formulation to robustly 
% extract the 4th-order tangent column by column.
% ============================================================

%% Number of Voigt components
if dim == 2
    nVoigt = 3;
elseif dim == 3
    nVoigt = 6;
else
    error('Unsupported dimension.');
end

%% Initialize Tangent
Cvoigt = zeros(nVoigt,nVoigt);

%% Principal strain trace and its derivative
trEps = sum(lambdaEig);

if strcmpi(branch,'plus')
    dtr = 0.5 * (1 + trEps/sqrt(trEps^2 + eta^2));
elseif strcmpi(branch,'minus')
    dtr = 0.5 * (1 - trEps/sqrt(trEps^2 + eta^2));
else
    error('Unknown spectral branch.');
end

%% ============================================================
% Construct tangent column by column (Kinematic Perturbation)
% ============================================================
for j = 1:nVoigt
    %% --------------------------------------------------------
    % Engineering strain perturbation (Basis Tensor H)
    % ---------------------------------------------------------
    H = voigtBasisTensor(j,dim);
    
    %% --------------------------------------------------------
    % Fréchet Derivative of Spectral Strain Function
    % ---------------------------------------------------------
    % Computes the exact change in the active/passive strain tensor 
    % (Df) resulting from the unit basis strain H.
    Df = zeros(3,3);
    
    for a = 1:3
        na = V(:,a);
        Pa = na*na'; % Eigen-projection tensor a
        
        for b = 1:3
            nb = V(:,b);
            Pb = nb*nb'; % Eigen-projection tensor b
            
            %% -----------------------------------------------
            % Divided Difference Evaluation
            % -----------------------------------------------
            if a == b
                % Diagonal components: No rotation coupling.
                % The derivative is simply the scalar derivative of the 
                % smoothed Macaulay bracket f'(lambda_a).
                coeff = spectralDerivative(lambdaEig(a),eta,branch);
            else
                % Off-diagonal components: Eigenvector rotation.
                % Represents how the principal axes shear against each other.
                deltaLambda = lambdaEig(a) - lambdaEig(b);
                scale = max([1; abs(lambdaEig(a)); abs(lambdaEig(b))]);
                
                % Check for distinct eigenvalues
                if abs(deltaLambda) > 1e-10 * scale
                    % Standard divided difference
                    coeff = (fEig(a) - fEig(b)) / deltaLambda;
                else
                    % ------------------------------------------------
                    % Repeated Eigenvalue Limit (L'Hôpital's Rule)
                    % ------------------------------------------------
                    % IMPORTANT: If the material is in a state of isotropic 
                    % tension/compression (lambda_a = lambda_b), the divided 
                    % difference yields 0/0 (a singularity). 
                    % By L'Hôpital's rule, the limit of the divided difference 
                    % as lambda_a -> lambda_b is simply the scalar derivative 
                    % f'(lambda). This unconditionally stabilizes the tangent 
                    % under hydrostatic stress states.
                    lambdaAvg = 0.5*(lambdaEig(a) + lambdaEig(b));
                    coeff = spectralDerivative(lambdaAvg,eta,branch);
                end
            end
            
            %% -----------------------------------------------
            % Accumulate Fréchet derivative contribution
            % -----------------------------------------------
            % Df = Sum_{a,b} coeff_ab * P_a * H * P_b
            Df = Df + coeff * Pa * H * Pb;
        end
    end
    
    %% --------------------------------------------------------
    % Derivative of active/passive stress
    % ---------------------------------------------------------
    % d(sigma_pm) = lambda * d(tr_pm)/d(tr) * tr(H) * I + 2 * mu * Df
    dsigma = lambdaLam*dtr*trace(H)*eye(3) + 2*mu*Df;
    
    %% --------------------------------------------------------
    % Convert to engineering stress notation (Column Mapping)
    % ---------------------------------------------------------
    Cvoigt(:,j) = stressTensorToVoigt(dsigma,dim);
end
end