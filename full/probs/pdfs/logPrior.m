function [log_prior, d_thetaR] = logPrior(thetaR, thetaP)
    %logPrior: This function computes the logarithm of the prior distribution
    %          for the system parameters. The prior consists of a mixture of Beta and 
    %          Gamma distributions depending on whether the parameter is bounded or 
    %          unbounded. If requested, the function also returns the derivative of the 
    %          log-prior with respect to the random parameters in thetaR.
    %
    % Inputs:
    %   thetaR   - Array of random system parameters for which the prior is evaluated
    %   thetaP   - Cell array containing the prior parameters (alphas, betas, ub, lb)
    %              for the Beta and Gamma distributions. ub and lb represent the 
    %              upper and lower bounds of each parameter, respectively.
    %
    % Outputs:
    %   log_prior - The log-prior for the given parameters
    %   d_thetaR  - The derivative of the log-prior with respect to thetaR
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % If thetaR is a cell, concatenate it into a single array
    if iscell(thetaR)
        thetaR = [thetaR{:}];
    end

    % Extract the prior parameters (alphas, betas, upper and lower bounds)
    [alphas, betas, ub, lb] = deal(thetaP{:});
    alphas = [alphas{:}];
    betas = [betas{:}];
    ub = [ub{:}];
    lb = [lb{:}];

    % Compute the width for bounded parameters
    width = ub - lb;

    % Identify unbounded and bounded parameters
    gs = isinf(ub);  % Unbounded (Gamma)
    bs = ~gs;        % Bounded (Beta)

    % Normalize the bounded parameters and shift the unbounded ones
    thetaR(:, bs) = (thetaR(:, bs) - lb(bs)) ./ width(bs);
    thetaR(:, gs) = thetaR(:, gs) - lb(gs);

    % If derivatives are requested
    if nargout > 1
        % Compute log-Beta and log-Gamma values and their derivatives
        [log_beta, dbetas]  = logBeta(thetaR(:, bs), alphas(bs), betas(bs));
        [log_gamma, dgammas] = logGamma(thetaR(:, gs), alphas(gs), betas(gs));

        % Adjust the derivatives for the bounded parameters
        dbetas = dbetas ./ width(bs);

        % Initialize the derivative array
        d_thetaR = zeros(size(thetaR));

        % Assign the derivatives for bounded and unbounded parameters
        d_thetaR(:, bs) = dbetas;
        d_thetaR(:, gs) = dgammas;
    else
        % Compute the log-Beta and log-Gamma values without derivatives
        log_beta  = logBeta(thetaR(:, bs), alphas(bs), betas(bs));
        log_gamma = logGamma(thetaR(:, gs), alphas(gs), betas(gs));
    end

    % Adjust log-Beta values for the bounded parameters
    log_beta = log_beta - log(width(bs));

    % Initialize the log-prior array
    log_prior = zeros(size(thetaR));

    % Assign the log-prior values for bounded and unbounded parameters
    log_prior(:, bs) = log_beta;
    log_prior(:, gs) = log_gamma;

    % Sum the log-prior values across dimensions
    log_prior = sum(log_prior, 2);
end