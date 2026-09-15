% ==========================================================
% Phase-Field Fracture Solver: PF-COMPAS
% AT2 Staggered
% Mesh-based formulation
%
% Preprocessing / Initialization / Adaptive Load Loop
%
% Author: Pranjal Saxena, IIT Kanpur
% Units: N, mm
% ==========================================================
clear;
clc;
tic;

%% =========================================================
% 1. MODEL DEFINITION
% ==========================================================
% The AT2 model relies heavily on the length scale parameter l0, 
% which governs the width of the diffuse crack. Gc is the critical 
% energy release rate (fracture toughness). xk is a small stabilization 
% parameter to prevent numerical singularity in fully damaged elements.
material.E  = 210000;       % Young's modulus [N/mm^2]
material.nu = 0.30;         % Poisson ratio
material.Gc = 2.7;          % Fracture energy [N/mm]
material.l0 = 0.04;         % Length scale [mm]
material.xk = 1e-8;         % Residual stiffness kappa

% ----------------------------------------------------------
% Phase-field split (Tension-Compression Asymmetry)
% 'none'     : Isotropic (computationally fast, physically limited)
% 'amor'     : Volumetric-Deviatoric decomposition
% 'spectral' : Principal strain decomposition (Miehe)
% ----------------------------------------------------------
model.split = 'amor';

% ----------------------------------------------------------
% Degradation Function g(phi)
% Couples the phase-field variable to the elastic strain energy.
% ----------------------------------------------------------
model.degradation = 'quadratic';

% =========================================================
% 2. LOAD MESH
% ==========================================================
meshFile = 'mesh2D_SENT_w_02_h_02.mat';
mesh = loadMesh(meshFile);

% Construct the assembly pattern for optimized sparse matrix allocation
assemblyPattern = buildAssemblyPattern(mesh);

nnodes = mesh.nodes.n;
nelem  = mesh.elem.n;
ndofU   = mesh.meta.dim * nnodes;
ndofPhi = nnodes;
model.dim = mesh.meta.dim;

debug.plotMesh = false;
if debug.plotMesh
    plot_mesh(mesh.nodes.coord, mesh.elem.conn, 'Q4', '-.');
end

% =========================================================
% 3. INITIALIZE CRACK AND BOUNDARY CONDITIONS 
% ==========================================================
geometry.tolerance = 1e-8;
x = mesh.nodes.x;
y = mesh.nodes.y;
tolGeom = geometry.tolerance;

% ----------------------------------------------------------
% Identify the initial SENT crack (y=0, x<0)
% ----------------------------------------------------------
crack_nodes1 = find(abs(y-0.01) < tolGeom & x < tolGeom);
crack_nodes2 = find(abs(y+0.01) < tolGeom & x < tolGeom);
crack_nodes = unique([crack_nodes1;crack_nodes2]);

% ---------------------------------------------------------
% Identify geometric boundaries for Dirichlet constraints
% ----------------------------------------------------------
top_nodes = find(abs(y - max(y)) < tolGeom);
bottom_nodes = find(abs(y - min(y)) < tolGeom);

% Pin one node in x-direction to prevent rigid-body translation
fix_node = find(abs(x - max(x)) < tolGeom & abs(y - min(y)) < tolGeom);

% ---------------------------------------------------------
% Displacement DOFs Mapping
% ----------------------------------------------------------
if model.dim == 2
    bc.u.topY = 2*top_nodes;
    bc.u.bottomY = 2*bottom_nodes;
    bc.u.fixX = 2*fix_node - 1;
elseif model.dim == 3
    bc.u.topY = 3*top_nodes - 1;
    bc.u.bottomY = 3*bottom_nodes - 1;
    bc.u.fixX = 3*fix_node - 2;
else
    error('Unsupported spatial dimension: %d', model.dim);
end

bc.u.loadedDofs = bc.u.topY(:);
bc.u.constrainedDofs = unique([bc.u.topY(:); bc.u.bottomY(:); bc.u.fixX(:)]);
bc.u.freeDofs = setdiff((1:ndofU)', bc.u.constrainedDofs);

% =========================================================
% PHASE-FIELD BOUNDARY CONDITIONS
% ==========================================================
% Enforce phi=1 (fully broken) at the initial crack nodes
bc.phi.fixedDofs = crack_nodes(:);
bc.phi.freeDofs = setdiff((1:ndofPhi)', bc.phi.fixedDofs);

% =========================================================
% 5. LOAD CONTROL (Displacement-driven)
% ==========================================================
% lambda is the dimensionless load parameter (0 to 1).
load.maxDisplacement = 0.01;     % [mm]
load.lambdaStart    = 0.0;
load.lambdaEnd      = 1.0;
load.dLambdaInitial = 0.02;
load.dLambdaMin     = 1e-8;
load.dLambdaMax     = 0.02;

% =========================================================
% 6. QUADRATURE
% ==========================================================
quad.type      = 'GAUSS';
quad.order     = 3;
quad.edgeOrder = 3;

% Standard element and edge quadrature points for numerical integration
[quad.W, quad.Q] = quadrature(quad.order, quad.type, model.dim);
[quad.edgeW, quad.edgeQ] = Gauss1D(quad.edgeOrder);

if model.dim == 3
    quad.faceOrder = 3;
    [quad.faceW,quad.faceQ] = quadrature(quad.faceOrder, quad.type, 2);
    quad.nFaceGP = numel(quad.faceW);
end
quad.nGP     = numel(quad.W);
quad.nEdgeGP = numel(quad.edgeW);

% =========================================================
% 7. SOLVER / NUMERICAL PARAMETERS
% ==========================================================
% Strict staggered limits to ensure accurate energy minimization
solver.maxStaggeredIter = 30;
solver.maxNewtonU   = 20;
solver.maxNewtonPhi = 20;

% Tolerances for internal Newton loops
solver.tolNewtonU      = 1e-9;
solver.tolNewtonPhi    = 1e-9;
solver.tolIncrementU   = 1e-9;
solver.tolIncrementPhi = 1e-9;

% Staggered loop convergence tolerances
solver.tolStaggeredU   = 1e-2;
solver.tolStaggeredPhi = 1e-6;

% Phase-field bounds (to maintain 0 <= phi <= 1)
solver.phiBoundTolerance = 1e-12;

% Targets for the adaptive load-stepping heuristic
solver.targetStaggeredIter = 15;
solver.targetNewtonU       = 10;
solver.targetNewtonPhi     = 10;
solver.targetDeltaPhi       = 0.02;
solver.maxAcceptedDeltaPhi  = 0.85;

solver.verbose = true;

% =========================================================
% 8. PARALLEL COMPUTING SETTINGS
% ==========================================================
parallel.useParpool = true;
parallel.useParfor  = true;
parallel.nWorkers = 4;
parallel.profile = 'local';
parallel.limitToAvailableWorkers = true;

if parallel.useParpool
    hasParallelToolbox = license('test', 'Distrib_Computing_Toolbox');
    if ~hasParallelToolbox
        warning('Parallel Computing Toolbox not available. Running in serial mode.');
        parallel.useParpool = false;
    else
        try
            localProfile = parallel.cluster.Local('NumWorkers', parallel.nWorkers);
            maxLocalWorkers = localProfile.NumWorkers;
        catch
            maxLocalWorkers = parallel.nWorkers;
        end
        if parallel.limitToAvailableWorkers
            parallel.nWorkers = min(parallel.nWorkers, maxLocalWorkers);
        end
        
        pool = gcp('nocreate');
        if isempty(pool)
            fprintf('\nStarting MATLAB parallel pool...\n');
            pool = parpool(parallel.profile, parallel.nWorkers);
        else
            fprintf('\nExisting parallel pool detected.\n');
            if pool.NumWorkers ~= parallel.nWorkers
                warning('Existing pool has %d workers. Existing pool will be used.', pool.NumWorkers);
            end
        end
        parallel.nWorkers = pool.NumWorkers;
        fprintf('Parallel processing ENABLED. Workers: %d\n\n', parallel.nWorkers);
    end
else
    fprintf('\nParallel processing DISABLED.\n\n');
end

% ===============================================
% 9. INITIALIZE PHASE FIELD FROM DISCRETE CRACK
% ===============================================
phi0 = initializePhaseField(mesh, crack_nodes, material, quad, parallel, assemblyPattern);

% =========================================================
% 10. SOLUTION STATE ALLOCATION
% ==========================================================
u   = zeros(ndofU,1);
phi = phi0;
H   = zeros(nelem,quad.nGP); % Strain history field for irreversibility

% [Output/ParaView Initializations omitted for brevity in comments - logic remains identical]
output.storeU          = true;
output.storePhi        = true;
output.storeStress     = true;
output.storeStressGP   = false;
output.storeHistory    = false;
output.solutionEvery = 10;
output.writeParaView = true;
output.saveResults = true;
output.baseName = 'SENT_AmorQuad_w_02_h_02_lo_04';
timeStamp = char(datetime('now','Format','yyyy_MM_dd_HH_mm_ss'));
output.resultsFolder = fullfile([output.baseName '_' timeStamp]);
if ~exist(output.resultsFolder,'dir'), mkdir(output.resultsFolder); end
output.vtuFolder = fullfile(output.resultsFolder,'paraview');
if output.writeParaView && ~exist(output.vtuFolder,'dir'), mkdir(output.vtuFolder); end

if output.writeParaView
    pvdFile = fullfile(output.vtuFolder, [output.baseName '.pvd']);
    initializePVDCollection(pvdFile);
else
    pvdFile = '';
end
if output.writeParaView
    pointData = struct();
    pointData.phi = phi;
    U3 = reshapeDisplacementToNodal(u, model.dim, nnodes);
    pointData.displacement = U3;
    pointData.Umag = sqrt(sum(U3.^2,2));
    cellData = struct();
    cellData.sigma_xx = zeros(nelem,1); cellData.sigma_yy = zeros(nelem,1);
    cellData.sigma_xy = zeros(nelem,1); cellData.vonMises = zeros(nelem,1);
    cellData.Eta_u = zeros(nelem,1); cellData.Eta_phi = zeros(nelem,1);
    cellData.Eta_mix = zeros(nelem,1);
    vtuFile = fullfile(output.vtuFolder, sprintf('%s_%04d.vtu', output.baseName,0));
    writeVTUStep(vtuFile, mesh, pointData, cellData);
    appendPVDEntry(pvdFile, vtuFile, load.lambdaStart);
end

resultsCapacity = 100;
nAcceptedMax    = resultsCapacity;
nStored = ceil(resultsCapacity/output.solutionEvery) + 1;

results = struct();
results.meta.problem   = 'AT2_phase_field';
results.meta.elemType = mesh.meta.elemType;
results.meta.dim      = model.dim;
results.meta.meshFile = meshFile;
results.meta.date     = datetime("now");
results.meta.material = material;
results.meta.model    = model;

results.time.acceptedStep = []; results.time.trialStep = [];
results.time.lambda = []; results.time.deltaLambda = []; results.time.load = [];
results.force.RF = []; results.energy.elastic = []; results.energy.fracture = []; results.energy.total = [];
results.indicator.u = []; results.indicator.phi = []; results.indicator.mix = [];
results.energyBalance.deltaElastic = []; results.energyBalance.deltaFracture = [];
results.energyBalance.deltaInternal = []; results.energyBalance.deltaExternal = [];
results.energyBalance.residual = []; results.energyBalance.absError = []; results.energyBalance.relativeError = [];

results.solver.nStaggered = []; results.solver.nNewtonU = []; results.solver.nNewtonPhi = [];
results.solution.step = []; results.solution.lambda = []; results.solution.load = [];
if output.storeU, results.solution.u = {}; end
if output.storePhi, results.solution.phi = {}; end
if output.storeStress, results.stress.sigma_avg = {}; results.stress.vonMises_avg = {}; end

lambda = load.lambdaStart; lambdaEnd = load.lambdaEnd; dLambda = load.dLambdaInitial;
step = 0; acceptedIndex = 1; storedIndex = 1; rejectedSteps = 0;

acceptedStep = 0; resultIndex = 1;
results.time.step(1) = 0; results.time.lambda(1) = lambda; results.time.load(1) = 0;
results.solution.step(1) = 0; results.solution.lambda(1) = lambda; results.solution.load(1) = 0;
if output.storeU, results.solution.u{1} = u; end
if output.storePhi, results.solution.phi{1} = phi; end

u0 = u; phi0 = phi; H0 = H;

[Eelastic0,Efracture0,Etotal0] = postprocessAT2Energy(mesh, u0, phi0, material, model, quad, parallel);
RF0 = computeReactionForce(mesh, u0, phi0, model, material, quad, 'top', 2);

results.time.step(acceptedIndex) = 0; results.time.lambda(acceptedIndex) = lambda;
results.time.load(acceptedIndex) = lambda * load.maxDisplacement;
results.force.RF(acceptedIndex) = RF0;
results.energy.elastic(acceptedIndex) = Eelastic0;
results.energy.fracture(acceptedIndex) = Efracture0;
results.energy.total(acceptedIndex) = Etotal0;

results.energyBalance.deltaElastic(acceptedIndex)  = 0; results.energyBalance.deltaFracture(acceptedIndex) = 0;
results.energyBalance.deltaInternal(acceptedIndex) = 0; results.energyBalance.deltaExternal(acceptedIndex) = 0;
results.energyBalance.residual(acceptedIndex) = 0; results.energyBalance.absError(acceptedIndex) = 0;
results.energyBalance.relativeError(acceptedIndex) = 0;

uAccepted = u; phiAccepted = phi; HAccepted = H;

%% =========================================================
% ADAPTIVE LOAD LOOP
% ==========================================================
while lambda < lambdaEnd - 1e-14
    %% PROPOSE NEXT LOAD STEP
    dLambdaTrial = min(dLambda,lambdaEnd-lambda);
    lambdaTrial = lambda + dLambdaTrial;
    step = step + 1;
    
    fprintf('\n============================================================\n');
    fprintf(' TRIAL LOAD STEP %d\n current lambda = %.8e\n trial lambda   = %.8e\n', step, lambda, lambdaTrial);
    fprintf('============================================================\n');

    % Backups to permit step rejection
    uTrialOld = uAccepted; phiTrialOld = phiAccepted; HTrialOld = HAccepted;

    % Update Boundary Conditions
    bcTrial = updateDisplacementBC(bc, lambdaTrial, load);
    
    % THERMODYNAMIC IRREVERSIBILITY
    % Prevent crack healing (phi_dot >= 0) by setting the minimum allowed phi
    % in the current step to the accepted phi from the previous step.
    bcTrial.phiMin = phiTrialOld;

    %% STAGGERED SOLVE
    [uTrial,phiTrial,HTrial,infoStag] = solveStaggeredLoadStep(mesh, uTrialOld, phiTrialOld, HTrialOld, bcTrial, material, model, quad, parallel, assemblyPattern, solver);

    %% CHECK TRIAL CONVERGENCE
    if ~infoStag.converged
        fprintf('\n  *** TRIAL STEP FAILED ***\n');
        rejectedSteps = rejectedSteps + 1;
        dLambda = 0.5*dLambda; % Cut step size in half
        
        if dLambda < load.dLambdaMin
            finalizePVDCollection(pvdFile);    
            resultsFile = fullfile(output.resultsFolder, [output.baseName '_results.mat']);
            save(resultsFile, 'results', '-v7.3');
            error('Adaptive load stepping failed. delta lambda reached minimum.');
        end
        clear uTrial phiTrial HTrial infoStag
        continue;
    end

    %% POSTPROCESS ENERGIES & FORCES
    [Eelastic,Efracture,Etotal] = postprocessAT2Energy(mesh, uTrial, phiTrial, material, model, quad, parallel);
    RF = computeReactionForce(mesh, uTrial, phiTrial, model, material, quad, 'top', 2);

    %% INCREMENTAL ENERGY BALANCE
    % Evaluates the accuracy of the staggered scheme by ensuring the external 
    % work done equals the change in internal elastic and fracture energies.
    uDispOld = lambda * load.maxDisplacement;
    uDispTrial = lambdaTrial * load.maxDisplacement;
    deltaU = uDispTrial - uDispOld;
    RFold = results.force.RF(acceptedIndex);
    
    deltaEelastic = Eelastic - results.energy.elastic(acceptedIndex);
    deltaEfracture = Efracture - results.energy.fracture(acceptedIndex);
    deltaEinternal = deltaEelastic + deltaEfracture;
    deltaWexternal = 0.5*(RFold + RF)*deltaU;
    
    energyResidual = deltaEinternal - deltaWexternal;
    energyAbsError = abs(energyResidual);
    energyScale = max([abs(deltaEinternal), abs(deltaWexternal), 1e-14]);
    energyRelativeError = energyAbsError/energyScale;
    
    %% SECONDARY CONVERGENCE CHECK (Excessive Damage Evolution)
    deltaPhi = phiTrial(bc.phi.freeDofs) - phiTrialOld(bc.phi.freeDofs);
    maxDeltaPhi = max(abs(deltaPhi));
   
    if maxDeltaPhi > solver.maxAcceptedDeltaPhi
        fprintf('  Trial rejected: excessive phase-field evolution\n');
        dLambda = 0.5*dLambdaTrial;
        if dLambda < load.dLambdaMin
            finalizePVDCollection(pvdFile);    
            resultsFile = fullfile(output.resultsFolder, [output.baseName '_results.mat']);
            save(resultsFile, 'results', '-v7.3');
            error('Adaptive load stepping failed. Excessive damage propagation.');
        end
        rejectedSteps = rejectedSteps + 1;
        clear uTrial phiTrial HTrial infoStag
        continue;
    end

    %% ACCEPT TRIAL STEP
    u = uTrial; phi = phiTrial; H = HTrial; lambda = lambdaTrial;
    uAccepted = u; phiAccepted = phi; HAccepted = H;

    %% ADAPTIVE STEP SIZE CONTROL (Heuristic)
    % Dynamically adjusts the load step based on how difficult the previous 
    % step was to solve, aiming to maintain a target number of Newton iterations.
    avgNewtonU = infoStag.totalNewtonU / max(infoStag.iterations,1);
    avgNewtonPhi = infoStag.totalNewtonPhi / max(infoStag.iterations,1);
    
    Dsolver = max([infoStag.iterations / solver.targetStaggeredIter, avgNewtonU / solver.targetNewtonU, avgNewtonPhi / solver.targetNewtonPhi]);
    Dphi = maxDeltaPhi / solver.targetDeltaPhi;
    difficulty = max(Dsolver,Dphi);
    
    factor = min(1.5,max(0.5, difficulty^(-0.5)));
    dLambda = min(load.dLambdaMax, max(load.dLambdaMin, dLambdaTrial*factor));

    %% EXPENSIVE POSTPROCESSING (Verfürth Indicators & Stresses)
    isFinalStep = abs(lambda - load.lambdaEnd) < 1e-14;
    storeSolution = mod(acceptedIndex,output.solutionEvery) == 0 || isFinalStep;
    needStress = output.storeStress && storeSolution || output.writeParaView;

    if needStress
        stress_avg = postprocessElementStress(mesh, u, phi, material, model, quad, parallel);
        vonMises = computeVonMises2D(stress_avg);
    else
        stress_avg = []; vonMises = [];
    end

    % Calculate rigorous Verfürth a posteriori error estimators
    [etaUGlobal,eta_elemU,~] = u2DResidualIndicator(mesh, u, phi, material, model, quad, dir_nodes);
    [etaPhiGlobal,eta_elemPhi,~] = phi2DResidualIndicator(mesh, phi, H, material, model, quad, crack_nodes);
    etaMixGlobal = sqrt(etaUGlobal.^2 + etaPhiGlobal.^2);
    eta_elemMix = sqrt(eta_elemU.^2 + eta_elemPhi.^2); 

    %% STORE ACCEPTED-STEP HISTORY
    acceptedIndex = acceptedIndex + 1;
    results.time.acceptedStep(acceptedIndex) = acceptedIndex;
    results.time.trialStep(acceptedIndex) = step;
    results.time.lambda(acceptedIndex) = lambda;
    results.time.deltaLambda(acceptedIndex) = dLambdaTrial;  
    results.force.RF(acceptedIndex) = RF;
    results.energy.elastic(acceptedIndex) = Eelastic;
    results.energy.fracture(acceptedIndex) = Efracture;
    results.energy.total(acceptedIndex) = Etotal;
    
    results.indicator.u(acceptedIndex) = etaUGlobal;
    results.indicator.phi(acceptedIndex) = etaPhiGlobal;
    results.indicator.mix(acceptedIndex) = etaMixGlobal;
    results.adapt.difficulty(acceptedIndex) = difficulty;
    
    results.energyBalance.deltaInternal(acceptedIndex) = deltaEinternal;
    results.energyBalance.deltaExternal(acceptedIndex) = deltaWexternal;
    results.energyBalance.relativeError(acceptedIndex) = energyRelativeError; 

    if storeSolution
        storedIndex = storedIndex + 1;
        results.solution.acceptedStep(storedIndex) = acceptedIndex;
        if output.storeU, results.solution.u{storedIndex} = u; end
        if output.storePhi, results.solution.phi{storedIndex} = phi; end
    end    

    %% PARAVIEW OUTPUT
    if output.writeParaView
        pointData.phi = phi;
        U3 = reshapeDisplacementToNodal(u, model.dim, nnodes);
        pointData.displacement = U3;
        pointData.Umag = sqrt(sum(U3.^2,2));
        if ~isempty(stress_avg)
            cellData.sigma_xx = stress_avg(:,1); cellData.sigma_yy = stress_avg(:,2);
            cellData.sigma_xy = stress_avg(:,3); cellData.vonMises = vonMises;
        end
        if ~isempty(eta_elemU)
            cellData.Eta_u = eta_elemU; cellData.Eta_phi = eta_elemPhi; cellData.Eta_mix = eta_elemMix;
        end
        vtuFile = fullfile(output.vtuFolder, sprintf('%s_%04d.vtu', output.baseName, acceptedIndex));
        writeVTUStep(vtuFile, mesh, pointData, cellData);
        appendPVDEntry(pvdFile, vtuFile, lambda*load.maxDisplacement);
    end    

    fprintf('\n  ACCEPTED STEP\n  lambda        = %.8e\n  staggered     = %d\n  energy error  = %.3e %%\n', lambda, infoStag.iterations, 100*energyRelativeError);
    
    clear stress_avg vonMises etaUGlobal etaPhiGlobal etaMixGlobal uTrial phiTrial HTrial
end
if output.writeParaView, finalizePVDCollection(pvdFile); end

%% =========================================================
% FINAL RESULTS SAVE & PLOTTING
% ==========================================================
if output.saveResults
    resultsFile = fullfile(output.resultsFolder, [output.baseName '_results.mat']);
    save(resultsFile, 'results', '-v7.3');
    
    scriptFile = matlab.desktop.editor.getActiveFilename;
    if ~isempty(scriptFile)
        [~, scriptName, scriptExt] = fileparts(scriptFile);
        scriptCopy = fullfile(output.resultsFolder, [output.baseName '_' scriptName scriptExt]);
        copyfile(scriptFile, scriptCopy);
    end
end

fig1 = figure;
uDisp = results.time.lambda.*load.maxDisplacement;
plot(uDisp, results.force.RF, 'LineWidth',2);
xlabel('\delta'); ylabel('Reaction force'); grid on; title('RF–\delta');
savefig(fig1, fullfile(output.resultsFolder, [output.baseName '_RF_delta.fig']));

fig2 = figure;
plot(uDisp, results.energy.elastic, 'b','LineWidth',2); hold on;
plot(uDisp, results.energy.fracture,'r','LineWidth',2);
plot(uDisp, results.energy.total, 'm--','LineWidth',2);
legend('Elastic','Fracture','Total'); xlabel('\delta'); ylabel('Energy'); grid on;
savefig(fig2, fullfile(output.resultsFolder, [output.baseName '_Energy.fig']));
toc;