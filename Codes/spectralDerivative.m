function df = spectralDerivative(lambda,eta,branch)
% ============================================================
% spectralDerivative - Derivative of Smoothed Principal Strains
%
% PHASE-FIELD HPC RELEVANCE:
% In the Miehe spectral split, the active and passive stresses are 
% functions of the positive and negative principal strains (lambda).
% To construct the exact consistent tangent moduli (C_tan), we must 
% differentiate these principal stresses with respect to the principal strains.
%
% Mathematically, this evaluates to:
%       d(sigma_plus_a) / d(lambda_a) = 2 * mu * f_plus'(lambda_a)
%
% If the standard absolute value function was used for the split, this 
% derivative would evaluate to the Heaviside step function (1 for lambda > 0, 
% 0 for lambda < 0), causing a discontinuous jump in the tangent matrix 
% that destroys Newton-Raphson convergence. 
%
% By differentiating our C^1-continuous smoothed Macaulay bracket instead, 
% this function returns a smooth, mathematically well-defined derivative 
% across the entire tension-compression transition zone.
%
% INPUT:
%   lambda : Principal strain eigenvalue
%   eta    : Dynamic smoothing parameter
%   branch : 'plus' (active) or 'minus' (passive)
%
% OUTPUT:
%   df     : Scalar derivative of the smoothed principal strain
% ============================================================

% Compute the hyperbolic smoothing root
root = sqrt(lambda^2 + eta^2);

if strcmpi(branch,'plus')
    % --------------------------------------------------------
    % Active (Tensile) Derivative
    % Limits: -> 1 as lambda >> 0
    %         -> 0 as lambda << 0
    %         -> 0.5 at lambda = 0
    % --------------------------------------------------------
    df = 0.5*(1 + lambda/root);
    
elseif strcmpi(branch,'minus')
    % --------------------------------------------------------
    % Passive (Compressive) Derivative
    % Limits: -> 0 as lambda >> 0
    %         -> 1 as lambda << 0
    %         -> 0.5 at lambda = 0
    % --------------------------------------------------------
    df = 0.5*(1 - lambda/root);
    
else
    error('Unknown spectral branch.');
end
end