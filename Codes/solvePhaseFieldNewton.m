function [phi,info] = ...
    solvePhaseFieldNewton( ...
    mesh,phi,H,bc, ...
    model,material,quad, ...
    parallel,assemblyPattern,solver)
% ============================================================
% solvePhaseFieldNewton - Active-Set Newton Solver for Damage
%
% PHASE-FIELD PHYSICS RELEVANCE:
% During this stage of the alternate minimization loop, the displacement 
% field (u) and the strain history field (H) are held strictly constant. 
% This renders the phase-field sub-problem strictly convex.
%
% THERMODYNAMIC IRREVERSIBILITY (ACTIVE-SET METHOD):
% While the history field H ensures the *driving force* never decreases, 
% numerical oscillations in the non-linear Newton steps can still 
% momentarily push phi below its previously accepted state (phiMin), 
% violating the thermodynamic condition that cracks cannot heal (phi_dot >= 0).
%
% Instead of using computationally stiff penalty methods, this framework 
% implements a highly elegant "Reduced-Space" or "Active-Set" Newton method. 
% It identifies nodes that are trying to heal (active DOFs) and explicitly 
% removes them from the iterative linear solve, guaranteeing exact, hard 
% enforcement of the irreversibility constraint.
% ============================================================

%% ============================================================
% INITIALIZE INFORMATION
% ============================================================
info = struct();
info.converged      = false;
info.iterations     = 0;
info.finalResidual  = inf;
info.finalIncrement = inf;
info.residualHistory  = zeros(solver.maxNewtonPhi,1);
info.incrementHistory = zeros(solver.maxNewtonPhi,1);

%% ============================================================
% FIXED CRACK DOFs (Geometric boundaries)
% ============================================================
fixedPhi = bc.phi.fixedDofs;
freePhi  = bc.phi.freeDofs;

%% ============================================================
% LOWER BOUND FOR IRREVERSIBILITY
% ============================================================
phiMin = bc.phiMin; % phi_old from the previously converged staggered step

if numel(phiMin) ~= numel(phi)
    error('bc.phiMin must have the same number of entries as phi.');
end

%% ============================================================
% INITIAL PHASE FIELD
% ============================================================
% Hard enforcement of irreversibility and initial macroscopic cracks 
% prior to beginning the Newton iterations.
phi = max(phi,phiMin);
phi(fixedPhi) = 1.0;

%% ============================================================
% INITIAL RESIDUAL SCALE
% ============================================================
R0 = [];

%% ============================================================
% NEWTON-RAPHSON LOOP
% ============================================================
for iter = 1:solver.maxNewtonPhi
    %% ========================================================
    % ASSEMBLE CURRENT STATE
    %
    % At iteration k: K_phi(phi^k) and R_phi(phi^k)
    % ========================================================
    [Kphi,Rphi] = ...
        assemblePhaseFieldSystem( ...
        mesh, ...
        phi, ...
        H, ...
        model,...
        material, ...
        quad, ...
        parallel, ...
        assemblyPattern);
        
    %% ========================================================
    % FREE DOF RESIDUAL
    % ========================================================
    Rfree = Rphi(freePhi);
    
    %% ========================================================
    % ACTIVE-SET DETECTION (The Irreversibility Check)
    %
    % Lower-bound constraint: phi >= phiMin
    %
    % If the Newton system is K * dphi = -R, then a positive residual 
    % (R > 0) implies that the solver will compute a negative increment 
    % (dphi < 0). 
    % 
    % If a node is already at its lower bound (phi <= phiMin + tol) AND 
    % R > 0, the unconstrained Newton step will force the crack to heal. 
    % This DOF is flagged as 'active' and will be locked.
    %
    % If R < 0, the solver wants to increase damage, which is physically 
    % valid, so the DOF remains 'inactive' (free to evolve).
    % ========================================================
    boundTol = solver.phiBoundTolerance;
    
    active = ...
        (phi(freePhi) <= phiMin(freePhi) + boundTol) ...
        & ...
        (Rfree > 0);
        
    %% ========================================================
    % PROJECTED RESIDUAL
    % ========================================================
    % The residual for actively healing DOFs is projected to exactly 0, 
    % artificially establishing "equilibrium" for these constrained nodes.
    Rprojected = Rfree;
    Rprojected(active) = 0;
    
    %% ========================================================
    % RESIDUAL NORM
    % ========================================================
    resNorm = norm(Rprojected);
    
    if iter == 1
        R0 = max(resNorm,1.0);
    end
    relRes = resNorm/R0;
    
    if iter == 1
        incNorm = inf;
    end
    
    info.residualHistory(iter) = relRes;
    
    if solver.verbose
        nActive = nnz(active);
        if isfinite(incNorm)
            fprintf('      Phi Newton %2d : R = %.3e, dPhi = %.3e, active = %d\n', iter, relRes, incNorm, nActive);
        else
            fprintf('      Phi Newton %2d : R = %.3e, active = %d\n', iter, relRes, nActive);
        end
    end
    
    %% ========================================================
    % DUAL-CRITERIA CONVERGENCE CHECK
    % ========================================================
    if iter > 1
        if relRes < solver.tolNewtonPhi && incNorm < solver.tolIncrementPhi
            info.converged      = true;
            info.iterations     = iter;
            info.finalResidual  = relRes;
            info.finalIncrement = incNorm;
            info.residualHistory  = info.residualHistory(1:iter);
            info.incrementHistory = info.incrementHistory(1:iter);
            return;
        end
    end
    
    %% ========================================================
    % REDUCED NEWTON SYSTEM (Inactive DOFs only)
    %
    %       K_ii dphi_i = -R_i
    %
    % Active DOFs receive a strict zero Newton increment.
    % ========================================================
    inactive = ~active;
    deltaFree = zeros(size(Rfree));
    
    if any(inactive)
        % Extract the reduced stiffness matrix and residual vector 
        % corresponding ONLY to nodes that are legally allowed to evolve.
        Kii = Kphi(freePhi(inactive), freePhi(inactive));
        Ri  = Rfree(inactive);
        
        % Solve the reduced system
        deltaFree(inactive) = -Kii \ Ri;
    end
    
    %% ========================================================
    % GLOBAL NEWTON INCREMENT
    % ========================================================
    dphi = zeros(size(phi));
    dphi(freePhi) = deltaFree;
    
    %% ========================================================
    % RAW UPDATE AND PROJECTION
    % ========================================================
    phiTrial = phi + dphi;
    
    % Enforce hard physical limits
    phiNew = max(phiTrial,phiMin);
    phiNew(fixedPhi) = 1.0;
    
    %% ========================================================
    % ACTUAL UPDATE
    %
    % IMPORTANT: We compute the increment norm based on the actual 
    % projected change (phiNew - phi) rather than the raw algorithmic 
    % step (dphi). This ensures mathematical consistency if the lower 
    % bound clipped a small oscillation.
    % ========================================================
    dphiActual = phiNew - phi;
    
    incNorm = norm(dphiActual(freePhi)) / max(norm(phiNew(freePhi)),1.0);
    info.incrementHistory(iter) = incNorm;
    
    % Accept the Newton iterate
    phi = phiNew;
end

%% ============================================================
% NEWTON FAILURE
% ============================================================
info.iterations = solver.maxNewtonPhi;
info.finalResidual = info.residualHistory(end);
info.finalIncrement = info.incrementHistory(end);
info.residualHistory = info.residualHistory(1:solver.maxNewtonPhi);
info.incrementHistory = info.incrementHistory(1:solver.maxNewtonPhi);

if solver.verbose
    fprintf('      Phi Newton FAILED: R = %.3e, dPhi = %.3e\n', info.finalResidual, info.finalIncrement);
end

end