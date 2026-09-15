function [u,info] = ...
    solveDisplacementNewton( ...
    mesh,u,phi, ...
    bc, ...
    material,model,quad, ...
    parallel,assemblyPattern,solver)
% ============================================================
% solveDisplacementNewton - Non-linear Mechanical Equilibrium
%
% PHASE-FIELD PHYSICS RELEVANCE:
% If the framework was strictly utilizing an isotropic (no-split) or 
% Amor (volumetric-deviatoric) energy split, the mechanical problem 
% with a frozen phi field would be purely linear, requiring only a 
% single matrix inversion. 
%
% However, the Miehe spectral split decomposes the strain tensor into 
% principal strains via eigenvalue decomposition. The principal directions 
% (eigenvectors) and magnitudes (eigenvalues) are highly non-linear 
% functions of the displacement field u. Therefore, the tangent stiffness 
% matrix K_u depends on the current displacement state u^k, necessitating 
% a robust Newton-Raphson iterative solver even when phi is held constant.
%
% Newton equation:
%       K_u(u^k) * du = -R_u(u^k)
% Update:
%       u^(k+1) = u^k + du
%
% IMPORTANT:
% The framework optimizes computational cost by avoiding a separate 
% "residual-only" assembly. The assembly at the beginning of the next 
% Newton iteration provides both K(u^(k+1)) and R(u^(k+1)).
% ============================================================

%% ============================================================
% INITIALIZE INFORMATION
% ============================================================
info = struct();
info.converged      = false;
info.iterations     = 0;
info.finalResidual  = inf;
info.finalIncrement = inf;
info.residualHistory  = zeros(solver.maxNewtonU,1);
info.incrementHistory = zeros(solver.maxNewtonU,1);

%% ============================================================
% PRESCRIBED DISPLACEMENT (Dirichlet BCs)
% ============================================================
% Ensure the current displacement vector strictly satisfies the 
% geometric boundary conditions for the current load step.
u(bc.u.constrainedDofs) = bc.u.constrainedValues;

%% ============================================================
% FREE DOFS
% ============================================================
free = bc.u.freeDofs;

%% ============================================================
% INITIAL RESIDUAL SCALE
% ============================================================
% We normalize the residual against the initial out-of-balance force (R0) 
% of the first iteration. max(R0, 1.0) prevents NaN divergence if the 
% system is already in perfect equilibrium (R0 ~ 0).
R0 = [];

%% ============================================================
% NEWTON-RAPHSON LOOP
% ============================================================
for iter = 1:solver.maxNewtonU
    
    %% ========================================================
    % ASSEMBLE CURRENT NEWTON STATE
    %
    % At iteration k:
    %       K = Tangent Stiffness Matrix K(u^k, phi)
    %       Rint = Internal Force Vector R(u^k, phi)
    % ========================================================
    [K,Rint] = ...
        assembleDisplacementSystem( ...
        mesh, ...
        u, ...
        phi, ...
        material, ...
        model, ...
        quad, ...
        parallel, ...
        assemblyPattern);
        
    %% ========================================================
    % CURRENT RESIDUAL
    % ========================================================
    Rfree = Rint(free);
    resNorm = norm(Rfree);
    
    %% --------------------------------------------------------
    % Reference residual normalization
    % --------------------------------------------------------
    if iter == 1
        R0 = max(resNorm,1.0);
    end
    relRes = resNorm/R0;
    
    %% ========================================================
    % INITIALIZE INCREMENT
    % ========================================================
    % Set to Inf on the first pass so that a zero initial residual 
    % requires at least one evaluated increment to declare true convergence.
    if iter == 1
        incNorm = inf;
    end
    
    %% ========================================================
    % STORE CURRENT RESIDUAL
    % ========================================================
    info.residualHistory(iter) = relRes;
    
    if solver.verbose
        if isfinite(incNorm)
            fprintf('      U Newton %2d : R = %.3e, dU = %.3e\n', iter, relRes, incNorm);
        else
            fprintf('      U Newton %2d : R = %.3e\n', iter, relRes);
        end
    end
    
    %% ========================================================
    % STRICT DUAL-CRITERIA CONVERGENCE CHECK
    %
    % The solver strictly demands that BOTH the normalized out-of-balance 
    % forces (relRes) AND the normalized displacement correction (incNorm) 
    % fall below their respective tolerances. 
    % This prevents "false convergence" where a highly stiff tangent matrix 
    % yields a tiny displacement update despite large residual forces.
    % ========================================================
    if iter > 1
        if relRes < solver.tolNewtonU && incNorm < solver.tolIncrementU
            info.converged      = true;
            info.iterations     = iter;
            info.finalResidual  = relRes;
            info.finalIncrement = incNorm;
            
            % Trim unused history pre-allocation
            info.residualHistory = info.residualHistory(1:iter);
            info.incrementHistory = info.incrementHistory(1:iter);
            return;
        end
    end
    
    %% ========================================================
    % SOLVE NEWTON SYSTEM (Reduced System)
    %
    %       K_ff * du_f = -R_f
    % ========================================================
    Kff = K(free,free);
    
    % Direct linear solve for the iterative correction
    duFree = -Kff \ Rfree;
    
    %% ========================================================
    % BUILD FULL INCREMENT & UPDATE
    % ========================================================
    du = zeros(size(u));
    du(free) = duFree;
    
    u = u + du; % u^(k+1) = u^k + du
    
    %% ========================================================
    % RE-ENFORCE PRESCRIBED VALUES
    %
    % Protects against floating-point corruption of constrained DOFs 
    % during the vector addition.
    % ========================================================
    u(bc.u.constrainedDofs) = bc.u.constrainedValues;
    
    %% ========================================================
    % INCREMENT CONVERGENCE MEASURE
    % ========================================================
    % Evaluates the L2 norm of the iterative correction relative to the 
    % total displacement magnitude.
    incNorm = norm(duFree) / max(norm(u(free)),1.0);
    info.incrementHistory(iter) = incNorm;
end

%% ============================================================
% NEWTON FAILURE
% ============================================================
info.iterations = solver.maxNewtonU;
info.finalResidual = info.residualHistory(end);
info.finalIncrement = info.incrementHistory(end);
info.residualHistory = info.residualHistory(1:solver.maxNewtonU);
info.incrementHistory = info.incrementHistory(1:solver.maxNewtonU);

if solver.verbose
    fprintf('      U Newton FAILED: R = %.3e, dU = %.3e\n', info.finalResidual, info.finalIncrement);
end

end