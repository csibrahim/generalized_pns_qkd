function [log_posterior, d_thetaE] = logPosterior(C, thetaA, thetaB, thetaE, thetaP)
    %logPosterior: Computes the logarithm of the posterior distribution, 
    %              which is the sum of the log-likelihood and the log-prior, 
    %              given the system parameters and observed counts. Optionally 
    %              calculates the gradient of the log-posterior with respect 
    %              to parameters in thetaE.
    %
    % Inputs:
    %   C       - Observed click counts
    %   thetaA  - Alice’s parameter set [lambdas, alpha, dAB]
    %   thetaB  - Bob’s parameter set [pa0, pa1, pc0, pc1, pd0, pd1, pe]
    %   thetaE  - Eve’s parameter set [dAE, pEB, k, Delta]
    %   thetaP  - Priors for Eve’s parameters [alphas, betas, ub, lb]
    %
    % Outputs:
    %   log_posterior - The log-posterior (sum of log-likelihood and log-prior)
    %   d_thetaE      - Derivative of the log-posterior with respect to 
    %                   parameters in thetaE
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    if nargout > 1
        % Compute log-likelihood and its derivatives
        [log_likelihood, d_likelihood] = logLikelihood(C, thetaA, thetaB, thetaE);

        % Compute log-prior and its derivatives
        [log_prior, d_prior] = logPrior(thetaE, thetaP);

        % Compute the derivative of the log-posterior
        d_thetaE = d_likelihood + d_prior;
    else
        % Compute log-likelihood without derivatives
        log_likelihood = logLikelihood(C, thetaA, thetaB, thetaE);

        % Compute log-prior without derivatives
        log_prior = logPrior(thetaE, thetaP);
    end

    % Compute the log-posterior as the sum of log-likelihood and log-prior
    log_posterior = log_likelihood + log_prior;
end
