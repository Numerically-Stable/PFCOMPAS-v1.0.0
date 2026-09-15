function [u,phi,H,info] = ...
    solveStaggeredLoadStep( ...
    mesh,uOld,phiOld,HOld, ...
    bc,material,model,quad, ...
    parallel,assemblyPattern,solver)
% ============================================================
% solveStaggeredLoadStep - Alternate Minimization Scheme
%
% PHASE-FIELD PHYSICS RELEVANCE:
% The total potential energy functional of the phase-field fracture problem 
% is highly non-convex with respect to the coupled variables (u, phi) 
% simultaneously. Attempting to solve them monolithically often leads to 
% severe Newton-Raphson divergence or convergence to non-physical local minima.
%
% However, the functional is strictly convex with respect to u (when phi is 
% fixed) and strictly convex with respect to phi (when u is fixed). 
% This framework utilizes an Alternate Minimization (Staggered) scheme, 
% passing back and forth between the displacement and damage fields until 
% a robust, global equilibrium is achieved within a single load increment.
%
% uOld   : previously accepted displacement
% phiOld : previously accepted phase field
% HOld   : previously accepted strain history field
% ============================================================

u   = uOld;
phi = phiOld;
H   = HOld;

info.converged = false;
info.iterations = 0;
info.totalNewtonU   = 0;
info.totalNewtonPhi = 0;
info.uChange   = inf;
info.phiChange = inf;

for iterStag = 1:solver.maxStaggeredIter
    if solver.verbose
        fprintf('\n');
        fprintf('   ----------------------------------------\n');
        fprintf('   STAGGERED ITERATION %d\n',iterStag);
        fprintf('   ----------------------------------------\n');
    end
    
    %% ========================================================
    % U SOLVE (Displacement Field)
    % ========================================================
    % Solve the balance of linear momentum holding the damage field 
    % (phi) completely constant. 
    uBefore = u;
    [uNew,infoU] = ...
        solveDisplacementNewton( ...
        mesh, ...
        u, ...
        phi, ...
        bc, ...
        material, ...
        model, ...
        quad, ...
        parallel, ...
        assemblyPattern, ...
        solver);
        
    info.UFinalResidual = infoU.finalResidual;
    
    % If the mechanical Newton solver fails, the entire staggered step fails.
    if ~infoU.converged
        info.converged = false;
        info.iterations = iterStag;
        return;
    end
    u = uNew;
    info.totalNewtonU = info.totalNewtonU + infoU.iterations;
    
    %% ========================================================
    % UPDATE HISTORY (Thermodynamic Irreversibility)
    % ========================================================
    % To enforce the thermodynamic condition that cracks cannot heal 
    % (damage must be monotonically increasing, phi_dot >= 0), we utilize 
    % a strain history field H. 
    % H = max( H_old, current_tensile_strain_energy )
    % The phase field is then driven by H rather than the instantaneous 
    % strain energy, preventing damage reduction during unloading phases.
    HNew = updateHistoryField( ...
        mesh, ...
        u, ...
        H, ...
        material, ...
        model, ...
        quad, ...
        parallel);
        
    %% ========================================================
    % PHI SOLVE (Phase Field)
    % ========================================================
    % Solve the damage evolution equation holding the mechanical 
    % displacement field (and thus HNew) completely constant.
    phiBefore = phi;
    [phiNew,infoPhi] = ...
        solvePhaseFieldNewton( ...
        mesh, ...
        phi, ...
        HNew, ...
        bc, ...
        model,...
        material, ...
        quad, ...
        parallel, ...
        assemblyPattern, ...
        solver);
        
    info.PhiFinalResidual = infoPhi.finalResidual;
    
    if ~infoPhi.converged
        info.converged = false;
        info.iterations = iterStag;
        return;
    end
    phi = phiNew;
    H = HNew;
    info.totalNewtonPhi = info.totalNewtonPhi + infoPhi.iterations;
    
    %% ========================================================
    % STAGGERED CHANGE (Convergence Metrics)
    % ========================================================
    % Calculate the relative L2 norm difference between consecutive 
    % staggered iterations to evaluate global equilibrium.
    duStag = ...
        norm(u-uBefore) / ...
        max(norm(u),1.0);
    dphiStag = ...
        norm(phi-phiBefore) / ...
        max(norm(phi),1.0);
        
    info.uChange   = duStag;
    info.phiChange = dphiStag;
    info.iterations = iterStag;
    
    if solver.verbose
        fprintf('   Staggered dU   = %.3e\n',duStag);
        fprintf('   Staggered dPhi = %.3e\n',dphiStag);
    end
    
    %% ========================================================
    % STAGGERED CONVERGENCE CHECK
    % ========================================================
    % In many practical phase-field setups, convergence of the displacement 
    % field (which governs the global force equilibrium) dominates. 
    % The threshold tolStaggeredU dictates when the alternate minimization 
    % has plateaued.
    if duStag < solver.tolStaggeredU 
        info.converged = true;
        if solver.verbose
            fprintf('   STAGGERED CONVERGED\n');
        end
        return;
    end
end

%% ============================================================
% STAGGERED FAILURE
% ============================================================
% If the maximum number of alternate minimizations is reached without 
% converging, the step is rejected, triggering the adaptive load 
% stepper in main.m to cut the increment size in half.
info.converged = false;
if solver.verbose
    fprintf( ...
        '   STAGGERED SOLVER FAILED after %d iterations.\n', ...
        solver.maxStaggeredIter);
end
end