% ==========================================================
% Phase-Field Fracture Solver: PF-COMPAS% ==========================================================
% Phase-Field Fracture Solver
% AT2 Staggered
% Mesh-based formulation
%
% Preprocessing / Initialization
%
% Author: Pranjal Saxena, IIT Kanpur
% Units: N, mm
% ==========================================================
clear;
clc;
tic;
%% Pre-processing %%
% =========================================================
% 1. MODEL DEFINITION
% ==========================================================
material.E  = 210000;       % Young's modulus [N/mm^2]
material.nu = 0.30;        % Poisson ratio
material.Gc = 2.7;         % Fracture energy [N/mm]
material.l0 = 0.1;       % Length scale [mm]
material.xk = 1e-8;        % Residual stiffness kappa
% ----------------------------------------------------------
% Phase-field split
%
%   'none'
%   'amor'
%   'spectral'
% ----------------------------------------------------------
model.split = 'amor';
% ----------------------------------------------------------
% Degradation
%
%   'quadratic'
%   'cubic'
%   'quartic'
%   'quintic'
%   'exponential'
%   'tanh'
%   'piecewiseCubic'
%   'piecewiseCubic2'
%   'augmentedTanh'
% ----------------------------------------------------------
model.degradation = 'quadratic';
% =========================================================
% 2. LOAD MESH
% ==========================================================
meshFile = 'mesh2D_SENT_w_02_h_02.mat';
mesh = loadMesh(meshFile);
assemblyPattern = buildAssemblyPattern(mesh);
nnodes = mesh.nodes.n;
nelem  = mesh.elem.n;
ndofU   = mesh.meta.dim * nnodes;
ndofPhi = nnodes;
% storing dimension into model
model.dim = mesh.meta.dim;
% Optional mesh plot
debug.plotMesh = false;
if debug.plotMesh
    plot_mesh( ...
        mesh.nodes.coord, ...
        mesh.elem.conn, ...
        'Q4', ...
        '-.');
end
% =========================================================
% 3. INITIALIZE CRACK AND BOUNDARY CONDITIONS 
% ==========================================================
geometry.tolerance = 1e-8;
x = mesh.nodes.x;
y = mesh.nodes.y;
tolGeom = geometry.tolerance;
% ----------------------------------------------------------
% Current SENT crack:
%
% y = 0
% x < 0
% ----------------------------------------------------------
crack_nodes1 = find( ...
    abs(y-0.01) < tolGeom & ...
    x < tolGeom);
crack_nodes2 = find( ...
    abs(y+0.01) < tolGeom & ...
    x < tolGeom);
crack_nodes = unique([crack_nodes1;crack_nodes2]);
%% ---------------------------------------------------------
% Identify boundary nodes for SENT
% ----------------------------------------------------------
top_nodes = find( ...
    abs(y - max(y)) < tolGeom);
bottom_nodes = find( ...
    abs(y - min(y)) < tolGeom);
% Pin one node in x-direction to remove rigid-body motion
fix_node = find( ...
    abs(x - max(x)) < tolGeom & ...
    abs(y - min(y)) < tolGeom);
dir_nodes = unique([ ...
    top_nodes(:);
    bottom_nodes(:);
    fix_node
]);
%% ---------------------------------------------------------
% Displacement DOFs
% ----------------------------------------------------------
% 2D:
%
%   x DOF = 2*i - 1
%   y DOF = 2*i
% 3D:
%   x DOF = 3*i - 2
%   y DOF = 3*i - 1
%   z DOF = 3*i
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
%% ---------------------------------------------------------
% Loaded displacement DOFs
% ----------------------------------------------------------
bc.u.loadedDofs = bc.u.topY(:);
%% ---------------------------------------------------------
% All constrained displacement DOFs
% ----------------------------------------------------------
bc.u.constrainedDofs =unique([bc.u.topY(:);bc.u.bottomY(:);bc.u.fixX(:)]);
%% ---------------------------------------------------------
% Free displacement DOFs
% ----------------------------------------------------------
bc.u.freeDofs =setdiff((1:ndofU)',bc.u.constrainedDofs);
%% =========================================================
% PHASE-FIELD BOUNDARY CONDITIONS
% ==========================================================
% Initial crack is prescribed as fully broken
bc.phi.fixedDofs = crack_nodes(:);
% Remaining phase-field DOFs
bc.phi.freeDofs = setdiff((1:ndofPhi)',bc.phi.fixedDofs);
% =========================================================
% 5. LOAD CONTROL
% ==========================================================
% Normalized load parameter:
%       lambda = 0  -> unloaded
%       lambda = 1  -> maximum displacement
% Physical displacement:
%       delta = lambda * maxDisplacement
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
% Element quadrature
[quad.W, quad.Q] = quadrature(quad.order,quad.type,model.dim);
% Edge quadrature
[quad.edgeW, quad.edgeQ] = Gauss1D(quad.edgeOrder);
% ----------------------------------------------------------
% 2D face quadrature for 3D
% This gives the tensor-product 3x3 rule 
% ----------------------------------------------------------
if model.dim == 3
    quad.faceOrder = 3;
    [quad.faceW,quad.faceQ] = quadrature(quad.faceOrder,quad.type,2);
    quad.nFaceGP = numel(quad.faceW);
end
quad.nGP     = numel(quad.W);
quad.nEdgeGP = numel(quad.edgeW);
% =========================================================
% 7. SOLVER / NUMERICAL PARAMETERS
% ==========================================================
solver.maxStaggeredIter = 30;
solver.maxNewtonU   = 20;
solver.maxNewtonPhi = 20;
% Newton residual tolerances
solver.tolNewtonU      = 1e-9;
solver.tolNewtonPhi    = 1e-9;
% Newton increment tolerances
solver.tolIncrementU   = 1e-9;
solver.tolIncrementPhi = 1e-9;
% Staggered convergence
solver.tolStaggeredU   = 1e-2;
solver.tolStaggeredPhi = 1e-6;
% Phase-field active-set tolerance
solver.phiBoundTolerance = 1e-12;
% Targets for adaptivity
solver.targetStaggeredIter = 20;
solver.targetNewtonU       = 15;
solver.targetNewtonPhi     = 15;
solver.targetDeltaPhi       = 0.002;
solver.maxAcceptedDeltaPhi  = 0.05;
% Output each info during solve
solver.verbose = true;
% =========================================================
% 8. PARALLEL COMPUTING SETTINGS
% ==========================================================
% IMPORTANT:
% Parallelization will eventually be used mainly for
% element-level calculations.
% ----------------------------------------------------------
parallel.useParpool = true;
parallel.useParfor  = true;
% Number of workers to use if parpool is enabled
parallel.nWorkers = 8;
% Local profile
parallel.profile = 'local';
% Automatically reduce requested worker count if more
% workers are requested than are physically available.
parallel.limitToAvailableWorkers = true;
% Start parallel pool
if parallel.useParpool
    % Check whether Parallel Computing Toolbox is available
    hasParallelToolbox = license( 'test','Distrib_Computing_Toolbox');
    if ~hasParallelToolbox
        warning(['Parallel Computing Toolbox not available. ','Running in serial mode.']);
        parallel.useParpool = false;
    else
        % Determine available local workers        
        try
            localProfile = parallel.cluster.Local('NumWorkers', parallel.nWorkers);
            maxLocalWorkers = localProfile.NumWorkers;
        catch
            % Fallback if cluster object cannot be queried
            maxLocalWorkers = parallel.nWorkers;
        end
        if parallel.limitToAvailableWorkers
            parallel.nWorkers = min(parallel.nWorkers,maxLocalWorkers);
        end
        % --------------------------------------------------
        % Check for an existing pool
        % --------------------------------------------------
        pool = gcp('nocreate');
        if isempty(pool)
            fprintf('\n');
            fprintf('Starting MATLAB parallel pool...\n');
            fprintf('Requested workers: %d\n',parallel.nWorkers);
            pool = parpool(parallel.profile,parallel.nWorkers);
        else
            fprintf('\n');
            fprintf('Existing parallel pool detected.\n');
            fprintf('Current workers: %d\n',pool.NumWorkers);
            if pool.NumWorkers ~= parallel.nWorkers
                warning(['Existing pool has %d workers, ', ...
                         'requested %d workers. ', ...
                         'Existing pool will be used.'], ...
                         pool.NumWorkers, ...
                         parallel.nWorkers);
            end
        end
        parallel.nWorkers = pool.NumWorkers;
        fprintf('Parallel processing ENABLED.\n');
        fprintf('Workers in use: %d\n',parallel.nWorkers);
        fprintf('\n');
    end
else
    fprintf('\n');
    fprintf('Parallel processing DISABLED.\n');
    fprintf('Running in serial mode.\n');
    fprintf('\n');
end
% ===============================================
% 9. INITIALIZE PHASE FIELD FROM DISCRETE CRACK
% ===============================================
phi0 = initializePhaseField( ...
    mesh, ...
    crack_nodes, ...
    material,...
    quad,...
    parallel,...
    assemblyPattern);
debug.plotphi0 = false;
if debug.plotphi0
    plotPhaseField(mesh, phi0);
end
% =========================================================
% 10. SOLUTION STATE
% ==========================================================
u   = zeros(ndofU,1);
phi = phi0;
H   = zeros(nelem,quad.nGP);
%% =========================================================
% 11. OUTPUT SETTINGS
% ==========================================================
output.storeU          = true;
output.storePhi        = true;
output.storeStress     = true;
output.storeStressGP   = false;
output.storeHistory    = false;
% Retain complete solution fields only every N accepted steps.
output.solutionEvery = 10;
% ParaView:
% Write every accepted/converged load step immediately.
output.writeParaView = true;
% Save results.mat only once after the complete analysis.
output.saveResults = true;
output.baseName = 'Sample_Results';
timeStamp = char(datetime('now','Format','yyyy_MM_dd_HH_mm_ss'));
output.resultsFolder =fullfile([output.baseName '_' timeStamp]);
if ~exist(output.resultsFolder,'dir')
    mkdir(output.resultsFolder);
end
output.vtuFolder = fullfile(output.resultsFolder,'paraview');
if output.writeParaView &&~exist(output.vtuFolder,'dir')
    mkdir(output.vtuFolder);
end
%% =========================================================
% INITIALIZE PARAVIEW
% ==========================================================
if output.writeParaView
    pvdFile = fullfile(output.vtuFolder,[output.baseName '.pvd']);
    initializePVDCollection(pvdFile);
else
    pvdFile = '';
end
if output.writeParaView
    pointData = struct();
    pointData.phi = phi;
    U3 = reshapeDisplacementToNodal(u,model.dim,nnodes);   
    pointData.displacement = U3;
    pointData.Umag = sqrt(sum(U3.^2,2));
    cellData = struct();
    cellData.sigma_xx = zeros(nelem,1);
    cellData.sigma_yy = zeros(nelem,1);
    cellData.sigma_xy = zeros(nelem,1);
    cellData.vonMises = zeros(nelem,1);
    cellData.Eta_u = zeros(nelem,1);
    cellData.Eta_phi = zeros(nelem,1);
    cellData.Eta_mix = zeros(nelem,1);
    vtuFile = fullfile(output.vtuFolder,sprintf('%s_%04d.vtu',output.baseName,0));
    writeVTUStep(vtuFile,mesh,pointData,cellData);
    appendPVDEntry(pvdFile,vtuFile,load.lambdaStart);
end
%% =========================================================
% OUTPUT CAPACITY
% ==========================================================
resultsCapacity = 100;
nAcceptedMax    = resultsCapacity;
% Number of solution states initially retained.
nStored = ceil(resultsCapacity/output.solutionEvery) + 1;
fprintf('\n');
fprintf('Maximum accepted steps : %d\n',nAcceptedMax);
fprintf('Stored full states     : %d\n',nStored);
fprintf('Solution retention     : every %d steps\n',output.solutionEvery);
fprintf('\n');
%% =========================================================
% RESULTS INITIALIZATION
% ==========================================================
results = struct();
results.meta.problem   = 'AT2_phase_field';
results.meta.elemType = mesh.meta.elemType;
results.meta.dim      = model.dim;
results.meta.meshFile = meshFile;
results.meta.date     = datetime("now");
results.meta.material = material;
results.meta.model    = model;
%% =========================================================
% ACCEPTED-STEP HISTORY
% ==========================================================
results.time.acceptedStep = [];
results.time.trialStep    = [];
results.time.lambda       = [];
results.time.deltaLambda  = [];
results.time.load         = [];
results.force.RF = [];
results.energy.elastic  = [];
results.energy.fracture = [];
results.energy.total    = [];
results.indicator.u     = [];
results.indicator.phi   = [];
results.indicator.mix   = [];
results.energyBalance.deltaElastic  = [];
results.energyBalance.deltaFracture = [];
results.energyBalance.deltaInternal = [];
results.energyBalance.deltaExternal = [];
results.energyBalance.residual      = [];
results.energyBalance.absError      = [];
results.energyBalance.relativeError = [];
%% =========================================================
% SOLVER DIAGNOSTICS
% ==========================================================
results.solver.nStaggered = [];
results.solver.nNewtonU   = [];
results.solver.nNewtonPhi = [];
%% =========================================================
% STORED FULL SOLUTIONS
% ==========================================================
results.solution.step   = [];
results.solution.lambda = [];
results.solution.load   = [];
if output.storeU
    results.solution.u = {};
end
if output.storePhi
    results.solution.phi = {};
end
%% =========================================================
% STORED STRESS
% ==========================================================
if output.storeStress
    results.stress.sigma_avg = {};
    results.stress.vonMises_avg = {};
end
%% =========================================================
% ADAPTIVE LOAD STEPPING
% ==========================================================
lambda = load.lambdaStart;
lambdaEnd = load.lambdaEnd;
dLambda = load.dLambdaInitial;
step = 0;
acceptedIndex = 1;
% Counts only states retained in RAM.
storedIndex = 1;
rejectedSteps = 0;
%% =========================================================
% INITIAL STATE
% ==========================================================
acceptedStep = 0;
resultIndex  = 1;
results.time.step(1)   = 0;
results.time.lambda(1) = lambda;
results.time.load(1)   = 0;
results.solution.step(1)   = 0;
results.solution.lambda(1) = lambda;
results.solution.load(1)   = 0;
if output.storeU
    results.solution.u{1} = u;
end
if output.storePhi
    results.solution.phi{1} = phi;
end
% ------------------------------------------------
% Initial state is lambda = 0
% ------------------------------------------------
u0   = u;
phi0 = phi;
H0   = H;
%% ================================================================
% POSTPROCESS INITIAL ENERGY
% ================================================================
[Eelastic0,Efracture0,Etotal0] =postprocessAT2Energy(mesh,u0,phi0,material,model,quad,parallel);
%% ================================================================
% INITIAL REACTION FORCE
% ================================================================
RF0 = computeReactionForce(mesh,u0,phi0,model,material,quad,'top',2);
%% ================================================================
% INITIAL HISTORY
% This is the reference state from which all increments are
% calculated.
%% ================================================================
results.time.step(acceptedIndex) = 0;
results.time.lambda(acceptedIndex) = lambda;
results.time.load(acceptedIndex) = lambda * load.maxDisplacement;
results.force.RF(acceptedIndex) = RF0;
results.energy.elastic(acceptedIndex) = Eelastic0;
results.energy.fracture(acceptedIndex) = Efracture0;
results.energy.total(acceptedIndex) = Etotal0;
% ------------------------------------------------
% No energy increment exists at the initial state
% ------------------------------------------------
results.energyBalance.deltaElastic(acceptedIndex)  = 0;
results.energyBalance.deltaFracture(acceptedIndex) = 0;
results.energyBalance.deltaInternal(acceptedIndex) = 0;
results.energyBalance.deltaExternal(acceptedIndex) = 0;
results.energyBalance.residual(acceptedIndex) = 0;
results.energyBalance.absError(acceptedIndex) = 0;
results.energyBalance.relativeError(acceptedIndex) = 0;
%% =========================================================
% INITIAL ACCEPTED STATE FOR ADAPTIVE CONTINUATION
% ==========================================================
uAccepted   = u;
phiAccepted = phi;
HAccepted   = H;
%% =========================================================
% ADAPTIVE LOAD LOOP
% ==========================================================
while lambda < lambdaEnd - 1e-14
    %% ======================================================
    % PROPOSE NEXT LOAD STEP
    %% ======================================================
    dLambdaTrial = min(dLambda,lambdaEnd-lambda);
    lambdaTrial = lambda + dLambdaTrial;
    step = step + 1;
    fprintf('\n');
    fprintf('============================================================\n');
    fprintf(' TRIAL LOAD STEP %d\n',step);
    fprintf(' current lambda = %.8e\n',lambda);
    fprintf(' trial lambda   = %.8e\n',lambdaTrial);
    fprintf(' delta lambda   = %.8e\n',dLambdaTrial);
    fprintf(' prescribed displacement = %.8e mm\n',lambdaTrial * load.maxDisplacement);
    fprintf('============================================================\n');
    %% ======================================================
    % SAVE CURRENT ACCEPTED STATE
    % IMPORTANT:
    % These MUST NOT be modified if the trial step fails.
    %% ======================================================
    uTrialOld   = uAccepted;
    phiTrialOld = phiAccepted;
    HTrialOld   = HAccepted;
    %% ======================================================
    % UPDATE BC FOR TRIAL LOAD
    %% ======================================================
    bcTrial = updateDisplacementBC(bc,lambdaTrial,load);
    %% ======================================================
    % IRREVERSIBILITY
    %% ======================================================
    bcTrial.phiMin = phiTrialOld;
    %% ======================================================
    % STAGGERED SOLVE
    %% ======================================================
    [uTrial,phiTrial,HTrial,infoStag] = ...
        solveStaggeredLoadStep( ...
        mesh, ...
        uTrialOld, ...
        phiTrialOld, ...
        HTrialOld, ...
        bcTrial, ...
        material, ...
        model, ...
        quad, ...
        parallel, ...
        assemblyPattern, ...
        solver);
    %% ======================================================
    % CHECK TRIAL CONVERGENCE
    %% ======================================================
    if ~infoStag.converged
        fprintf('\n');
        fprintf('  *** TRIAL STEP FAILED ***\n');
        rejectedSteps = rejectedSteps + 1;
        dLambda = 0.5*dLambda;
        fprintf('  Reducing delta lambda to %.4e\n',dLambda);
        if dLambda < load.dLambdaMin         
            finalizePVDCollection(pvdFile);    
            resultsFile = fullfile(output.resultsFolder, ...
                [output.baseName '_results.mat']);      
            fprintf('\n');
            fprintf('Saving final results...\n');
            fprintf('  %s\n',resultsFile);       
            save(resultsFile,'results','-v7.3');
            fprintf('Final results saved successfully.\n');
            error(['Adaptive load stepping failed.\n' ...
                 'delta lambda reached minimum %.3e\n' ...
                 'Current lambda = %.8e'],load.dLambdaMin,lambda);
        end
        % IMPORTANT:
        % Do NOT modify u_old, phi_old or H_old.
        % Retry from the last accepted state.
        clear uTrial phiTrial HTrial infoStag
        continue;
    end
    %% ======================================================
    % POSTPROCESS ACCEPTED TRIAL STATE
    %% ======================================================
    [Eelastic,Efracture,Etotal] = ...
        postprocessAT2Energy( ...
        mesh, ...
        uTrial, ...
        phiTrial, ...
        material, ...
        model, ...
        quad, ...
        parallel);
    %% ======================================================
    % REACTION FORCE
    %% ======================================================
    RF = computeReactionForce( ...
        mesh, ...
        uTrial, ...
        phiTrial, ...
        model, ...
        material, ...
        quad, ...
        'top', ...
        2);
    %% ======================================================
    % INCREMENTAL ENERGY BALANCE
    %% ======================================================
    uDispOld = lambda * load.maxDisplacement;
    uDispTrial = lambdaTrial * load.maxDisplacement;
    deltaU = uDispTrial - uDispOld;
    RFold = results.force.RF(acceptedIndex);
    EelasticOld = results.energy.elastic(acceptedIndex);
    EfractureOld = results.energy.fracture(acceptedIndex);
    deltaEelastic = Eelastic - EelasticOld;
    deltaEfracture = Efracture - EfractureOld;
    deltaEinternal = deltaEelastic + deltaEfracture;
    deltaWexternal = 0.5*(RFold + RF)*deltaU;
    energyResidual = deltaEinternal - deltaWexternal;
    energyAbsError = abs(energyResidual);
    energyScale = max([abs(deltaEinternal),abs(deltaWexternal),1e-14]);
    energyRelativeError = energyAbsError/energyScale;    
    %% ================================================
    % Second level of check for convergence
    %% ================================================
    deltaPhi =phiTrial(bc.phi.freeDofs) - phiTrialOld(bc.phi.freeDofs);
    maxDeltaPhi = max(abs(deltaPhi));
    if maxDeltaPhi > solver.maxAcceptedDeltaPhi
        fprintf('  Trial rejected: excessive phase-field evolution\n');    
        dLambda = 0.5*dLambdaTrial;
        if dLambda < load.dLambdaMin
            finalizePVDCollection(pvdFile);    
            resultsFile = fullfile(output.resultsFolder, ...
                [output.baseName '_results.mat']);      
            fprintf('\n');
            fprintf('Saving final results...\n');
            fprintf('  %s\n',resultsFile);        
            save(resultsFile,'results','-v7.3');      
            fprintf('Final results saved successfully.\n');
            error(['Adaptive load stepping failed.\n' ...
                 'delta lambda reached minimum %.3e\n' ...
                 'Current lambda = %.8e'],load.dLambdaMin,lambda);
        end
        rejectedSteps = rejectedSteps + 1;    
        clear uTrial phiTrial HTrial infoStag    
        continue;
    end
    %% ======================================================
    % ACCEPT TRIAL STEP
    %% ======================================================
    u   = uTrial;
    phi = phiTrial;
    H   = HTrial;
    lambda = lambdaTrial;
    uAccepted   = u;
    phiAccepted = phi;
    HAccepted   = H;
    %% =========================================================
    % ADAPTIVE STEP SIZE CONTROL
    % ==========================================================        
    % ---------------------------------------------------------
    % Solver difficulty
    % ---------------------------------------------------------
    avgNewtonU = infoStag.totalNewtonU / max(infoStag.iterations,1);    
    avgNewtonPhi = infoStag.totalNewtonPhi / max(infoStag.iterations,1);   
    Dsolver = max([infoStag.iterations / solver.targetStaggeredIter, ...
        avgNewtonU / solver.targetNewtonU, ...
        avgNewtonPhi / solver.targetNewtonPhi]);  
    %% ---------------------------------------------------------
    % Phase-field difficulty
    % ---------------------------------------------------------    
    Dphi = maxDeltaPhi / solver.targetDeltaPhi;    
    %% ---------------------------------------------------------
    % Overall difficulty
    % ----------------------------------------------------------    
    difficulty = max(Dsolver,Dphi);
    %% ---------------------------------------------------------
    % Step-size factor
    % ---------------------------------------------------------    
    factor = difficulty^(-0.5);
    %% ---------------------------------------------------------
    % Limit adaptation
    % ---------------------------------------------------------    
    factor = min(1.5,max(0.5,factor));    
    %% ---------------------------------------------------------
    % Next step
    % ---------------------------------------------------------
    dLambda = min(load.dLambdaMax,max(load.dLambdaMin,dLambdaTrial*factor));
    %% ======================================================
    % DETERMINE WHETHER FULL STATE IS RETAINED
    % ======================================================
    isFinalStep = abs(lambda - load.lambdaEnd) < 1e-14;
    storeSolution = mod(acceptedIndex,output.solutionEvery) == 0 || isFinalStep;
    %% ======================================================
    % EXPENSIVE POSTPROCESSING
    %
    % Required if:
    %   1. storing the state
    %   2. ParaView output is requested
    %
    % Otherwise do not calculate it.
    %% ======================================================
    needStress = ...
        output.storeStress && storeSolution || ...
        output.writeParaView;
    %% ------------------------------------------------------
    % Stress
    %% ------------------------------------------------------
    if needStress
        stress_avg = ...
            postprocessElementStress( ...
            mesh, ...
            u, ...
            phi, ...
            material, ...
            model, ...
            quad, ...
            parallel);
        vonMises = ...
            computeVonMises2D(stress_avg);
    else
        stress_avg = [];
        vonMises = [];
    end

    %% ------------------------------------------------------
    % Error indicators
    %% ------------------------------------------------------
    [etaUGlobal,eta_elemU,~] = u2DResidualIndicator(mesh, ...
                                                    u, ...
                                                    phi, ...
                                                    material, ...
                                                    model, ...
                                                    quad,...
                                                    dir_nodes);
    [etaPhiGlobal,eta_elemPhi,~] = phi2DResidualIndicator(mesh, ...
                                                          phi, ...
                                                          H, ...
                                                          material, ...
                                                          model, ...
                                                          quad, ...
                                                          crack_nodes);
     etaMixGlobal = sqrt(etaUGlobal.^2 + etaPhiGlobal.^2);
     eta_elemMix = sqrt(eta_elemU.^2 + eta_elemPhi.^2); 
    %% ======================================================
    % STORE ACCEPTED-STEP HISTORY
    % These arrays contain EVERY converged load step.
    %% ======================================================
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
    results.adapt.nRejected = rejectedSteps;
    results.adapt.difficulty(acceptedIndex) = difficulty;    
    results.adapt.maxDeltaPhi(acceptedIndex) = maxDeltaPhi;
    results.adapt.nextDeltaLambda(acceptedIndex) = dLambda;
    %% ------------------------------------------------------
    % Incremental energy balance
    %% ------------------------------------------------------
    results.energyBalance.deltaElastic(acceptedIndex) = deltaEelastic;
    results.energyBalance.deltaFracture(acceptedIndex) = deltaEfracture;
    results.energyBalance.deltaInternal(acceptedIndex) = deltaEinternal;
    results.energyBalance.deltaExternal(acceptedIndex) = deltaWexternal;
    results.energyBalance.residual(acceptedIndex) = energyResidual;
    results.energyBalance.absError(acceptedIndex) = energyAbsError;
    results.energyBalance.relativeError(acceptedIndex) = energyRelativeError; 
    %% ------------------------------------------------------
    % Solver statistics
    %% ------------------------------------------------------
    results.solver.nStaggered(acceptedIndex) = infoStag.iterations;
    results.solver.nNewtonU(acceptedIndex) = infoStag.totalNewtonU;
    results.solver.nNewtonPhi(acceptedIndex) = infoStag.totalNewtonPhi;
    results.solver.finalUResidual(acceptedIndex) = infoStag.UFinalResidual;
    results.solver.finalPhiResidual(acceptedIndex) = infoStag.PhiFinalResidual;
    results.solver.uChange(acceptedIndex) = infoStag.uChange;
    results.solver.phiChange(acceptedIndex) = infoStag.phiChange; 
    %% ======================================================
    % STORE FULL SOLUTION
    % Independent of acceptedIndex.
    %% ======================================================
    if storeSolution
        storedIndex = storedIndex + 1;
        results.solution.acceptedStep(storedIndex) = acceptedIndex;
        results.solution.trialStep(storedIndex) = step;
        results.solution.lambda(storedIndex) = lambda;
        results.solution.load(storedIndex) = lambda*load.maxDisplacement ;
        %% --------------------------------------------------
        % Displacement
        %% --------------------------------------------------
        if output.storeU
            results.solution.u{storedIndex} = u;
        end
        %% --------------------------------------------------
        % Phase field
        %% --------------------------------------------------
        if output.storePhi
            results.solution.phi{storedIndex} = phi;
        end
        %% --------------------------------------------------
        % Stress
        %% --------------------------------------------------
        if output.storeStress
            results.stress.sigma_avg{storedIndex} = stress_avg;
            results.stress.vonMises_avg{storedIndex} = vonMises;
        end
    end    
    %% ======================================================
    % PARAVIEW OUTPUT
    % Every accepted load step.
    %% ======================================================
    if output.writeParaView
        %% --------------------------------------------------
        % Point data
        %% --------------------------------------------------
        pointData = struct();
        pointData.phi = phi;
        U3 = reshapeDisplacementToNodal( ...
            u, ...
            model.dim, ...
            nnodes);        
        pointData.displacement = U3;
        pointData.Umag = sqrt(sum(U3.^2,2));
        %% --------------------------------------------------
        % Cell data
        %% --------------------------------------------------
        cellData = struct();
        if ~isempty(stress_avg)
            cellData.sigma_xx = stress_avg(:,1);
            cellData.sigma_yy = stress_avg(:,2);
            cellData.sigma_xy = stress_avg(:,3);
            cellData.vonMises = vonMises;
        end
        if ~isempty(eta_elemU)
            cellData.Eta_u = eta_elemU;
            cellData.Eta_phi = eta_elemPhi;
            cellData.Eta_mix = eta_elemMix;
        end
        %% --------------------------------------------------
        % VTU filename
        %% --------------------------------------------------
        vtuFile = fullfile(output.vtuFolder,sprintf('%s_%04d.vtu', ...
            output.baseName,acceptedIndex));
        %% --------------------------------------------------
        % Write VTU
        %% --------------------------------------------------
        writeVTUStep(vtuFile,mesh,pointData,cellData);
        %% --------------------------------------------------
        % Update PVD
        %% --------------------------------------------------
        appendPVDEntry(pvdFile,vtuFile,lambda*load.maxDisplacement);
        clear pointData cellData U3
    end    
    %% ======================================================
    % PRINT ACCEPTED STEP
    %% ======================================================
    fprintf('\n');
    fprintf('  ACCEPTED STEP\n');
    fprintf('  lambda        = %.8e\n',lambda);
    fprintf('  delta lambda  = %.8e\n',dLambdaTrial);
    fprintf('  next delta    = %.8e\n',dLambda);
    fprintf('  staggered     = %d\n',infoStag.iterations);
    fprintf('  max DeltaPhi  = %.3e\n',maxDeltaPhi);
    fprintf('  energy error  = %.3e %%\n',100*energyRelativeError);
    %% ======================================================
    % RELEASE STEP-LOCAL DATA
    %% ======================================================
    clear stress_avg vonMises
    clear etaUGlobal etaPhiGlobal etaMixGlobal
    clear uTrial phiTrial HTrial
    clear Eelastic Efracture Etotal RF
    clear u_new phi_new H_new infoStag   
end
if output.writeParaView
    finalizePVDCollection(pvdFile);    
end
%% =========================================================
% FINALIZE RESULTS
% ==========================================================
fprintf('\n');
fprintf('============================================================\n');
fprintf(' ANALYSIS COMPLETE\n');
fprintf('============================================================\n');
%% =========================================================
% FINAL RESULTS SAVE
%% ==========================================================
if output.saveResults
    resultsFile = fullfile(output.resultsFolder,[output.baseName '_results.mat']);
    scriptFile = mfilename('fullpath');
    [~, scriptName, scriptExt] = fileparts(scriptFile);
    scriptCopy = fullfile(output.resultsFolder,[output.baseName '_' scriptName scriptExt]);
    fprintf('\n');
    fprintf('Saving final results...\n');
    fprintf('  %s\n',resultsFile);
    save(resultsFile,'results','-v7.3');
    fprintf('Final results saved successfully.\n');
    fprintf('Saving script copy...\n');    
    scriptFile = matlab.desktop.editor.getActiveFilename;    
    if isempty(scriptFile)
        warning('Could not determine the active MATLAB script filename.');
    else    
        [~, scriptName, scriptExt] = fileparts(scriptFile);   
        scriptCopy = fullfile(output.resultsFolder, ...
            [output.baseName '_' scriptName scriptExt]); 
        fprintf('  %s\n', scriptCopy);   
        [status, msg] = copyfile(scriptFile, scriptCopy);  
        if status
            fprintf('Script copy saved successfully.\n');
        else
            warning('Could not save script copy:\n%s', msg);
        end   
    end
end

%% ========================================================================
% BASIC QUANTITIES
% ========================================================================
% Load/displacement history
uDisp = results.time.lambda .* load.maxDisplacement;
% Reaction force
RF = results.force.RF;
% Energies
Eelastic  = results.energy.elastic;
Efracture = results.energy.fracture;
Etotal    = results.energy.total;
%% ========================================================================
% 1. REACTION FORCE - DISPLACEMENT
% ========================================================================
fig1 = figure;
plot(uDisp, RF, 'LineWidth', 2);
xlabel('\delta');
ylabel('Reaction force');
title('Reaction force-displacement response');
grid on;
box on;
adaptiveName = [output.baseName '_RF_delta'];
savefig(gcf,fullfile(output.resultsFolder,[adaptiveName '.fig']));
saveas(gcf,fullfile(output.resultsFolder,[adaptiveName '.svg']));

%% ========================================================================
% 2. ENERGY EVOLUTION
% ========================================================================
fig2 = figure;
plot(uDisp, Eelastic,'LineWidth', 2);
hold on;
plot(uDisp, Efracture,'LineWidth', 2);
plot(uDisp, Etotal,'--','LineWidth', 2);
legend('Elastic','Fracture','Total','Location','best');
xlabel('\delta');
ylabel('Energy');
title('Energy evolution');
grid on;
box on;
adaptiveName = [output.baseName '_Energy'];
savefig(gcf,fullfile(output.resultsFolder,[adaptiveName '.fig']));
saveas(gcf,fullfile(output.resultsFolder,[adaptiveName '.svg']));

%% ========================================================================
% 3. VERFÜRTH ERROR INDICATORS
% ========================================================================
fig3 = figure;
plot(uDisp,results.indicator.u,'LineWidth', 2);
hold on;
plot(uDisp,results.indicator.phi,'LineWidth', 2);
legend('\eta_u','\eta_\phi','Location','best');
xlabel('\delta');
ylabel('Error indicator');
title('Verfürth error indicators');
grid on;
box on;
adaptiveName = [output.baseName '_ErrorIndicators'];
savefig(gcf,fullfile(output.resultsFolder,[adaptiveName '.fig']));
saveas(gcf,fullfile(output.resultsFolder,[adaptiveName '.svg']));

%% ========================================================================
% 4. VERFÜRTH ERROR INDICATORS - LOG SCALE
% ========================================================================
fig4 = figure;
semilogy(uDisp,results.indicator.u,'LineWidth', 2);
hold on;
semilogy(uDisp,results.indicator.phi,'LineWidth', 2);
legend('\eta_u','\eta_\phi','Location','best');
xlabel('\delta');
ylabel('Error indicator');
title('Verfürth error indicators - logarithmic scale');
grid on;
box on;
adaptiveName = [output.baseName '_ErrorIndicators_Log'];
savefig(gcf,fullfile(output.resultsFolder,[adaptiveName '.fig']));
saveas(gcf,fullfile(output.resultsFolder,[adaptiveName '.svg']));

%% ========================================================================
% 5. ADAPTIVE LOAD STEP SIZE
% ========================================================================
fig5 = figure;
semilogy(uDisp,results.time.deltaLambda,'LineWidth', 2);
xlabel('\delta');
ylabel('\Delta\lambda');
title('Adaptive load increment');
grid on;
box on;
adaptiveName = [output.baseName '_DeltaLambda'];
savefig(gcf,fullfile(output.resultsFolder,[adaptiveName '.fig']));
saveas(gcf,fullfile(output.resultsFolder,[adaptiveName '.svg']));

%% ========================================================================
% 6. LOAD STEP SIZE + MAXIMUM PHASE-FIELD CHANGE
% ========================================================================
fig6 = figure;
yyaxis left
semilogy(uDisp,results.time.deltaLambda,'LineWidth', 2);
ylabel('\Delta\lambda');
yyaxis right
plot(uDisp,results.adapt.maxDeltaPhi,'LineWidth', 2);
ylabel('max |\Delta\phi|');
xlabel('\delta');
title('Adaptive load stepping and phase-field increment');
grid on;
box on;
adaptiveName = [output.baseName '_AdaptiveStep_PhiChange'];
savefig(gcf,fullfile(output.resultsFolder,[adaptiveName '.fig']));
saveas(gcf,fullfile(output.resultsFolder,[adaptiveName '.svg']));

%% ========================================================================
% 7. ADAPTIVE DIFFICULTY
% ========================================================================
fig7 = figure;
plot(uDisp,results.adapt.difficulty,'LineWidth', 2);
xlabel('\delta');
ylabel('Difficulty');
title('Adaptive stepping difficulty');
grid on;
box on;
adaptiveName = [output.baseName '_Difficulty'];
savefig(gcf,fullfile(output.resultsFolder,[adaptiveName '.fig']));
saveas(gcf,fullfile(output.resultsFolder,[adaptiveName '.svg']));

%% ========================================================================
% 8. INCREMENTAL ENERGY BALANCE
% ========================================================================
fig8 = figure;
plot(uDisp,results.energyBalance.deltaInternal,'LineWidth', 2);
hold on;
plot(uDisp,results.energyBalance.deltaExternal,'LineWidth', 2);
plot(uDisp,results.energyBalance.residual,'--','LineWidth', 2);
legend('\Delta E_{internal}','\Delta W_{external}','Energy residual','Location','best');
xlabel('\delta');
ylabel('Incremental energy');
title('Incremental energy balance');
grid on;
box on;
adaptiveName = [output.baseName '_EnergyBalance'];
savefig(gcf,fullfile(output.resultsFolder,[adaptiveName '.fig']));
saveas(gcf,fullfile(output.resultsFolder,[adaptiveName '.svg']));

%% ========================================================================
% 9. RELATIVE ENERGY BALANCE ERROR
% ========================================================================
fig9 = figure;
semilogy(uDisp,abs(results.energyBalance.relativeError),'LineWidth', 2);
xlabel('\delta');
ylabel('Relative energy error');
title('Relative energy balance error');
grid on;
box on;
adaptiveName = [output.baseName '_EnergyBalanceRelativeError'];
savefig(gcf,fullfile(output.resultsFolder,[adaptiveName '.fig']));
saveas(gcf,fullfile(output.resultsFolder,[adaptiveName '.svg']));

%% ========================================================================
% 10. STAGGERED ITERATIONS
% ========================================================================
fig10 = figure;
plot(uDisp,results.solver.nStaggered,'LineWidth', 2);
xlabel('\delta');
ylabel('Number of staggered iterations');
title('Staggered solver iterations');
grid on;
box on;
adaptiveName = [output.baseName '_StaggeredIterations'];
savefig(gcf,fullfile(output.resultsFolder,[adaptiveName '.fig']));
saveas(gcf,fullfile(output.resultsFolder,[adaptiveName '.svg']));

%% ========================================================================
% 11. NEWTON ITERATIONS
% ========================================================================
fig11 = figure;
plot(uDisp,results.solver.nNewtonU,'LineWidth', 2);
hold on;
plot(uDisp,results.solver.nNewtonPhi,'LineWidth', 2);
legend('Displacement','Phase field','Location','best');
xlabel('\delta');
ylabel('Number of Newton iterations');
title('Newton solver iterations');
grid on;
box on;
adaptiveName = [output.baseName '_NewtonIterations'];
savefig(gcf,fullfile(output.resultsFolder,[adaptiveName '.fig']));
saveas(gcf,fullfile(output.resultsFolder,[adaptiveName '.svg']));

%% ========================================================================
% 12. NONLINEAR SOLVER RESIDUALS
% ========================================================================
fig12 = figure;
semilogy(uDisp,abs(results.solver.finalUResidual),'LineWidth', 2);
hold on;
semilogy(uDisp,abs(results.solver.finalPhiResidual),'LineWidth', 2);
legend('Displacement residual','Phase-field residual','Location','best');
xlabel('\delta');
ylabel('Residual');
title('Final nonlinear residuals');
grid on;
box on;
adaptiveName = [output.baseName '_SolverResiduals'];
savefig(gcf,fullfile(output.resultsFolder,[adaptiveName '.fig']));
saveas(gcf,fullfile(output.resultsFolder,[adaptiveName '.svg']));

%% ========================================================================
% 13. SOLUTION CHANGES
% ========================================================================
fig13 = figure;
semilogy(uDisp,abs(results.solver.uChange),'LineWidth', 2);
hold on;
semilogy(uDisp,abs(results.solver.phiChange),'LineWidth', 2);
legend('||\Delta u||','||\Delta\phi||','Location','best');
xlabel('\delta');
ylabel('Solution change');
title('Solution increments');
grid on;
box on;
adaptiveName = [output.baseName '_SolutionChanges'];
savefig(gcf,fullfile(output.resultsFolder,[adaptiveName '.fig']));
saveas(gcf,fullfile(output.resultsFolder,[adaptiveName '.svg']));

%% ========================================================================
% 14. COMPLETE ANALYSIS SUMMARY
% ========================================================================
fig14 = figure;
tiledlayout(2,2);

% -------------------------------------------------------------------------
% Reaction force
% -------------------------------------------------------------------------
nexttile;
plot(uDisp,RF,'LineWidth', 2);
xlabel('\delta');
ylabel('Reaction force');
title('Load response');
grid on;
box on;

% -------------------------------------------------------------------------
% Energy
% -------------------------------------------------------------------------
nexttile;
plot(uDisp,Eelastic,'LineWidth', 2);
hold on;
plot(uDisp,Efracture,'LineWidth', 2);
plot(uDisp,Etotal,'--','LineWidth', 2);
xlabel('\delta');
ylabel('Energy');
title('Energy');
legend('Elastic','Fracture','Total','Location','best');
grid on;
box on;

% -------------------------------------------------------------------------
% Error indicators
% -------------------------------------------------------------------------
nexttile;
semilogy(uDisp,results.indicator.u,'LineWidth', 2);
hold on;
semilogy(uDisp,results.indicator.phi,'LineWidth', 2);
xlabel('\delta');
ylabel('Error indicator');
title('Verfürth indicators');
legend('\eta_u','\eta_\phi','Location','best');
grid on;
box on;

% -------------------------------------------------------------------------
% Energy balance
% -------------------------------------------------------------------------
nexttile;
semilogy(uDisp,abs(results.energyBalance.absError),'LineWidth', 2);
xlabel('\delta');
ylabel('Absolute error');
title('Energy balance');
grid on;
box on;

sgtitle('PFCOMPAS-v1.0.0 analysis summary');
adaptiveName = [output.baseName '_AnalysisSummary'];
savefig(gcf,fullfile(output.resultsFolder,[adaptiveName '.fig']));
saveas(gcf,fullfile(output.resultsFolder,[adaptiveName '.svg']));

%% ========================================================================
% FINISHED
% ========================================================================
fprintf('\n');
fprintf('============================================================\n');
fprintf('Analysis completed successfully.\n');
fprintf('============================================================\n');
fprintf('Results folder:\n%s\n',output.resultsFolder);
fprintf('============================================================\n');
fprintf('\n');
toc;
