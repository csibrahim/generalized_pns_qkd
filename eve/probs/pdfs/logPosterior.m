function [log_posterior, d_thetaE] = logPosterior(C, thetaA, thetaB, thetaE, thetaP)
    % logPosterior: This function computes the logarithm of the posterior 
    %               distribution, which is the sum of the log-likelihood and 
    %               the log-prior, given the system parameters and observed counts. 
    %               If requested, it also computes the derivative of the log-posterior 
    %               with respect to the parameters in thetaE.
    %
    % Inputs:
    %   C          - Array of observed counts for each outcome.
    %   thetaA     - Cell array of system parameters specific to the channel between 
    %                Alice and Bob, containing [lambdas, alpha, dAB].
    %   thetaB     - Cell array of system parameters for Bob's detectors, containing 
    %                [pa0, pa1, pc0, pc1, pd0, pd1, pe].
    %   thetaE     - Cell array of system parameters related to Eve’s eavesdropping, 
    %                containing [dAE, pEB, k, Delta].
    %   thetaP     - Cell array containing the prior parameters (alphas, betas, ub, lb).
    %
    % Outputs:
    %   log_posterior - The log-posterior, computed as the sum of the log-likelihood 
    %                   and log-prior.
    %   d_thetaE      - The gradient of the log-posterior with respect to the 
    %                   parameters in thetaE.
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    if (nargout > 1)
        % Compute log-likelihood and its gradient
        [log_likelihood, d_likelihood] = logLikelihood(C, thetaA, thetaB, thetaE);

        % Compute log-prior and its gradient
        [log_prior, d_prior] = logPrior(thetaE, thetaP);
    
        % Combine gradients for the log-posterior
        d_thetaE = d_likelihood + d_prior;
    else
        % Compute log-likelihood without gradient
        log_likelihood = logLikelihood(C, thetaA, thetaB, thetaE);

        % Compute log-prior without gradient
        log_prior = logPrior(thetaE, thetaP);
    end

    % Sum log-likelihood and log-prior to get log-posterior
    log_posterior = log_likelihood + log_prior;
end