function [g,dg,ddg] = degradationFunction(phi,model,material)
% ============================================================
% degradationFunction - Constitutive Stiffness Degradation
%
% PHASE-FIELD PHYSICS RELEVANCE:
% The degradation function g(phi) fundamentally couples the scalar damage 
% field to the mechanical displacement field. As phi transitions from 
% 0 (intact) to 1 (fully fractured), g(phi) monotonically decreases from 
% 1 to 0, physically degrading the material's load-bearing capacity.
%
% WHY SO MANY OPTIONS?
% The classical AT2 model uses a standard 'quadratic' degradation: (1-phi)^2. 
% However, this causes stiffness reduction to begin immediately upon loading, 
% eliminating the purely linear elastic phase and prematurely initiating damage. 
% Modern phase-field cohesive zone models (PF-CZM) utilize higher-order 
% (cubic, quartic) or asymptotic (tanh, exponential) functions to artificially 
% delay crack nucleation until a specific critical stress/strain threshold 
% is reached, bringing the model closer to true experimental behaviors.
%
% ALGORITHMIC RELEVANCE OF DERIVATIVES:
%   g   : Degrades the mechanical stiffness matrix (K_u) and stress tensor.
%   dg  : Acts as the energetic driving force for damage in the internal 
%         residual vector of the phase-field evolution equation (R_phi).
%   ddg : Strictly required for assembling the exact Jacobian (tangent 
%         stiffness matrix, K_phi) in the phase-field Newton-Raphson solver 
%         to guarantee quadratic convergence rates.
% ============================================================

%% ============================================================
% Residual stiffness (kappa)
% ============================================================
% As phi -> 1, g(phi) -> 0. If the structural stiffness matrix K_u were 
% allowed to reach strictly zero, the global mechanical system would become 
% singular and the linear solver would crash. Kappa (k ~ 10^-8) acts as an 
% artificial, infinitesimal residual stiffness to maintain a positive-definite 
% matrix across fully fractured domains.
if isfield(material,'kappa')
    kappa = material.kappa;
else
    kappa = 0;
end

%% ============================================================
% Check degradation choice
% ============================================================
if ~isfield(model,'degradation')
    error('material.degradation must be specified.');
end
name = lower(model.degradation);

%% ============================================================
% Degradation function definitions
% ============================================================
switch name
    %% ========================================================
    % 1. QUADRATIC (Standard AT2 Model)
    %
    % g = (1-phi)^2 + kappa
    % ========================================================
    case 'quadratic'
        g   = (1-phi)^2 + kappa;
        dg  = -2*(1-phi);
        ddg = 2;
        
    %% ========================================================
    % 2. CUBIC
    %
    % Flattens the slope of g(phi) near phi=0, delaying the onset 
    % of macroscopic softening.
    % ========================================================
    case 'cubic'
        g   = 3*(1-phi)^2 - 2*(1-phi)^3 + kappa;
        dg  = -6*(1-phi)  + 6*(1-phi)^2;
        ddg = 6 - 12*(1-phi);
        
    %% ========================================================
    % 3. QUARTIC
    % ========================================================
    case 'quartic'
        g   = 1 - 6*phi^2 + 8*phi^3 - 3*phi^4 + kappa;
        dg  = -12*phi + 24*phi^2 - 12*phi^3;
        ddg = -12 + 48*phi - 36*phi^2;
        
    %% ========================================================
    % 4. QUINTIC
    %
    % C2 smoothstep-type degradation. Ensures that both the first 
    % and second derivatives vanish smoothly at the boundaries, 
    % highly stabilizing the Newton solver.
    % ========================================================
    case 'quintic'
        g   = 1 - 10*phi^3 + 15*phi^4 - 6*phi^5 + kappa;
        dg  = -30*phi^2 + 60*phi^3 - 30*phi^4;
        ddg = -60*phi + 180*phi^2 - 120*phi^3;
        
    %% ========================================================
    % 5. EXPONENTIAL
    %
    % Creates a sharp, asymptotic drop in stiffness, closely mimicking 
    % brittle cleavage.
    % ========================================================
    case 'exponential'
        if isfield(model,'degradationParameter')
            a = model.degradationParameter;
        else
            a = 8;
        end
        expTerm = exp(-a*phi^2);
        g   = expTerm + kappa;
        dg  = -2*a*phi*expTerm;
        ddg = (-2*a + 4*a^2*phi^2)*expTerm;
        
    %% ========================================================
    % 6. TANH
    %
    % Smooth step transition shifted by xc.
    % ========================================================
    case 'tanh'
        if isfield(model,'degradationParameter')
            p = model.degradationParameter;
            xc = p(1);
            w  = p(2);
        else
            xc = 0.5;
            w  = 0.1;
        end
        z = (phi-xc)/w;
        sech2 = sech(z)^2;
        g   = 0.5*(1-tanh(z)) - 0.000045 + kappa;
        dg  = -0.5/w * sech2;
        ddg = 1/w^2 * sech2*tanh(z);
        
    %% ========================================================
    % 7. PIECEWISE CUBIC (B-Spline inspired)
    % ========================================================
    case 'piecewisecubic'
        if phi <= 0.01
            a = 10000; b = -200; c = 0; d = 1;
        elseif phi >= 0.99
            a = 10000; b = -29800; c = 29600; d = -9800;
        else
            a = 0; b = 0; c = -1; d = 1;
        end
        g   = a*phi^3 + b*phi^2 + c*phi + d + kappa;
        dg  = 3*a*phi^2 + 2*b*phi + c;
        ddg = 6*a*phi + 2*b;
        
    %% ========================================================
    % 8. PIECEWISE CUBIC 2
    % ========================================================
    case 'piecewisecubic2'
        if phi <= 0.2
            a = -7.5; b = -0.25; c = 0; d = 1;
        elseif phi >= 0.99
            a = 10000; b = -29800; c = 29600; d = -9800;
        else
            a = 0.5273; b = -0.9413; c = -0.6868; d = 1.1008;
        end
        g   = a*phi^3 + b*phi^2 + c*phi + d + kappa;
        dg  = 3*a*phi^2 + 2*b*phi + c;
        ddg = 6*a*phi + 2*b;
        
    %% ========================================================
    % 9. AUGMENTED TANH
    %
    % Highly customizable S-curve allowing for independent control 
    % of the nucleation threshold (m) and softening rate (w).
    % ========================================================
    case 'augmentedtanh'
        if isfield(model,'degradationParameter')
            p = model.degradationParameter;
            xc = p(1); w  = p(2); m  = p(3);
        else
            xc = 0.85; w  = 0.15; m  = 1.2;
        end
        %% ----------------------------------------------------
        % Warping
        % -----------------------------------------------------
        s   = 1 - (1-phi)^m;
        ds  = m*(1-phi)^(m-1);
        d2s = -m*(m-1)*(1-phi)^(m-2);
        
        %% ----------------------------------------------------
        % Normalized tanh
        % -----------------------------------------------------
        z     = (s-xc)/w;
        th    = tanh(z);
        sech2 = sech(z)^2;
        
        %% ----------------------------------------------------
        % Normalization constants
        % -----------------------------------------------------
        T0 = tanh((0-xc)/w);
        T1 = tanh((1-xc)/w);
        D  = T1-T0;
        
        %% ----------------------------------------------------
        % Degradation and exact analytical derivatives
        % -----------------------------------------------------
        g   = (T1-th)/D + kappa;
        dg  = -(sech2*ds)/(D*w);
        ddg = ( 2*ds^2*sech2*th/w^2 - d2s*sech2/w ) / D;
        
    %% ========================================================
    % UNKNOWN
    % ========================================================
    otherwise
        error('Unknown degradation function: %s', model.degradation);
end
end