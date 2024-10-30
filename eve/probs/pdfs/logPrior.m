function [log_prior, d_thetaE] = logPrior(thetaE, thetaP)
    %logPrior: This function computes the logarithm of the prior distribution 
    %          for the system parameters. The prior consists of a mixture of 
    %          Beta and Gamma distributions, depending on whether each parameter 
    %          is bounded or unbounded. If requested, it also computes the gradient 
    %          of the log-prior with respect to each parameter in thetaE.
    %
    % Inputs:
    %   thetaE - Array of system parameters related to Eve’s eavesdropping, for which 
    %            the prior is evaluated.
    %   thetaP - Cell array containing the prior parameters (alphas, betas, ub, lb)
    %            for the Beta and Gamma distributions. ub and lb represent the 
    %            upper and lower bounds for each parameter, respectively.
    %
    % Outputs:
    %   log_prior - The log-prior value for the given parameters.
    %   d_thetaE  - Cell array containing the derivative of the log-prior with 
    %               respect to each parameter in thetaE.
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % If thetaE is a cell array, concatenate it into a single array
    if iscell(thetaE)
        thetaE = [thetaE{:}];
    end

    % Extract the prior parameters (alphas, betas, upper and lower bounds)
    [alphas, betas, ub, lb] = deal(thetaP{:});
    width = ub - lb;  % Calculate the width for bounded parameters

    % Identify bounded (Beta) and unbounded (Gamma) parameters
    gs = isinf(ub);  % Unbounded parameters (Gamma)
    bs = ~gs;        % Bounded parameters (Beta)

    % Normalize the bounded parameters and shift the unbounded ones
    thetaE(:, bs) = (thetaE(:, bs) - lb(bs)) ./ width(bs);
    thetaE(:, gs) = thetaE(:, gs) - lb(gs);

    if nargout > 1
        % Compute log-Beta and log-Gamma values and their derivatives
        [log_beta, d_betas]   = logBeta(thetaE(:, bs), alphas(bs), betas(bs));
        [log_gamma, d_gammas] = logGamma(thetaE(:, gs), alphas(gs), betas(gs));
    
        % Adjust derivatives for the bounded parameters
        d_betas = d_betas ./ width(bs);
    
        % Initialize the gradient array
        d_thetaE = zeros(size(thetaE));
    
        % Assign the derivatives for bounded and unbounded parameters
        d_thetaE(:, bs) = d_betas;
        d_thetaE(:, gs) = d_gammas;
    else
        % If no derivatives are requested, only compute the log-probabilities
        log_beta  = logBeta(thetaE(:, bs), alphas(bs), betas(bs));
        log_gamma = logGamma(thetaE(:, gs), alphas(gs), betas(gs));

    end

    % Adjust log-Beta values for bounded parameters
    log_beta = log_beta - log(width(bs));

    % Initialize the log-prior array
    log_prior = zeros(size(thetaE));

    % Assign log-prior values for bounded and unbounded parameters
    log_prior(:, bs) = log_beta;
    log_prior(:, gs) = log_gamma;

    % Sum the log-prior values across dimensions
    log_prior = sum(log_prior, 2);
end