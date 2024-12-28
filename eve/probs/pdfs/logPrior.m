function [log_prior, d_thetaE] = logPrior(thetaE, thetaP)
    %logPrior: Computes the logarithm of the prior distribution for system 
    %          parameters. The prior is a mixture of Beta distributions for
    %          bounded parameters and Gamma distributions for semi-bounded 
    %          parameters. Optionally calculates the gradient of the log-prior 
    %          with respect to parameters in thetaE.
    %
    % Inputs:
    %   thetaE - Eve’s parameter set [dAE, pEB, k, Delta]
    %   thetaP - Priors for Eve’s parameters [alphas, betas, ub, lb]
    %
    % Outputs:
    %   log_prior - The log-prior value for the given parameters
    %   d_thetaE  - The derivative of the log-prior with respect to thetaE
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Flatten thetaE if provided as a cell array
    if iscell(thetaE)
        thetaE = [thetaE{:}];
    end

    % Extract prior parameters
    [alphas, betas, ub, lb] = deal(thetaP{:});

    % Compute the width for bounded parameters
    width = ub - lb;

    % Identify semi-bounded and bounded parameters
    gs = isinf(ub);  % Semi-bounded (Gamma) parameters
    bs = ~gs;        % Bounded (Beta) parameters

    % Normalize bounded parameters and shift semi-bounded ones
    thetaE(:, bs) = (thetaE(:, bs) - lb(bs)) ./ width(bs);
    thetaE(:, gs) = thetaE(:, gs) - lb(gs);

    % If derivatives are requested
    if nargout > 1
        % Compute log-probabilities and their derivatives
        [log_beta, d_betas]   = logBeta(thetaE(:, bs), alphas(bs), betas(bs));
        [log_gamma, d_gammas] = logGamma(thetaE(:, gs), alphas(gs), betas(gs));

        % Adjust derivatives for normalized bounded parameters
        d_betas = d_betas ./ width(bs);

        % Initialize the derivative array
        d_thetaE = zeros(size(thetaE));

        % Assign derivatives for bounded and semi-bounded parameters
        d_thetaE(:, bs) = d_betas;
        d_thetaE(:, gs) = d_gammas;
    else
        % Compute only log-probabilities if derivatives are not requested
        log_beta  = logBeta(thetaE(:, bs), alphas(bs), betas(bs));
        log_gamma = logGamma(thetaE(:, gs), alphas(gs), betas(gs));
    end

    % Adjust log-probabilities for bounded parameters
    log_beta = log_beta - log(width(bs));

    % Initialize log-prior array
    log_prior = zeros(size(thetaE));

    % Assign log-prior values for bounded and semi-bounded parameters
    log_prior(:, bs) = log_beta;
    log_prior(:, gs) = log_gamma;

    % Sum the log-prior values across dimensions
    log_prior = sum(log_prior, 2);
end
