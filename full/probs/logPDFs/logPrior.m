function [log_prior, d_thetaR] = logPrior(thetaR, thetaP)
    %logPrior: Computes the logarithm of the prior distribution for system 
    %          parameters. The prior is a mixture of Beta distributions for
    %          bounded parameters and Gamma distributions for semi-bounded 
    %          parameters. Optionally calculates the gradient of the log-prior 
    %          with respect to parameters in thetaE.
    %
    % Inputs:
    %   thetaR   - Random variables' values for which the prior is evaluated
    %   thetaP   - Prior parameters of random variables [alphas, betas, ub, lb] 
    %
    % Outputs:
    %   log_prior - The log-prior value for the given parameters
    %   d_thetaR  - The derivative of the log-prior with respect to thetaR
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Flatten thetaR if provided as a cell array
    if iscell(thetaR)
        thetaR = [thetaR{:}];
    end

    % Extract prior parameters
    [alphas, betas, ub, lb] = deal(thetaP{:});

    % Flatten prior parameters
    alphas = [alphas{:}];
    betas = [betas{:}];
    ub = [ub{:}];
    lb = [lb{:}];

    % Compute the width for bounded parameters
    width = ub - lb;

    % Identify semi-bounded and bounded parameters
    gs = isinf(ub);  % Semi-bounded (Gamma)
    bs = ~gs;        % Bounded (Beta)

    % Normalize bounded parameters and shift semi-bounded ones
    thetaR(:, bs) = (thetaR(:, bs) - lb(bs)) ./ width(bs);
    thetaR(:, gs) = thetaR(:, gs) - lb(gs);

    % If derivatives are requested
    if nargout > 1
        % Compute log-probabilities and their derivatives
        [log_beta, dbetas]  = logBeta(thetaR(:, bs), alphas(bs), betas(bs));
        [log_gamma, dgammas] = logGamma(thetaR(:, gs), alphas(gs), betas(gs));

        % Adjust derivatives for normalized bounded parameters
        dbetas = dbetas ./ width(bs);

        % Initialize the derivative array
        d_thetaR = zeros(size(thetaR));

        % Assign derivatives for bounded and semi-bounded parameters
        d_thetaR(:, bs) = dbetas;
        d_thetaR(:, gs) = dgammas;
    else
        % Compute only log-probabilities if derivatives are not requested
        log_beta  = logBeta(thetaR(:, bs), alphas(bs), betas(bs));
        log_gamma = logGamma(thetaR(:, gs), alphas(gs), betas(gs));
    end

    % Adjust log-probabilities for bounded parameters
    log_beta = log_beta - log(width(bs));

    % Initialize log-prior array
    log_prior = zeros(size(thetaR));

    % Assign log-prior values for bounded and semi-bounded parameters
    log_prior(:, bs) = log_beta;
    log_prior(:, gs) = log_gamma;

    % Sum the log-prior values across dimensions
    log_prior = sum(log_prior, 2);
end