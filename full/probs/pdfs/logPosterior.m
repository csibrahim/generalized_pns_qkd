function [log_posterior, d_thetaR] = logPosterior(thetaR, thetaF, varR, varF, C, thetaP)
    %logPosterior: This function computes the logarithm of the posterior 
    %              distribution, which is the sum of the log-likelihood and 
    %              the log-prior, given the system parameters and observed counts. 
    %              If requested, it also computes the derivative of the log-posterior 
    %              with respect to the random parameters in thetaR.
    %
    % Inputs:
    %   thetaR        - Array of random system parameters
    %   thetaF        - Array of fixed system parameters
    %   varR          - Cell array of variable names corresponding to the random parameters
    %   varF          - Cell array of variable names corresponding to the fixed parameters
    %   C             - Array of observed counts for each outcome
    %   thetaP        - Cell array containing the prior parameters (alphas, betas, ub, lb)
    %
    % Outputs:
    %   log_posterior - The log-posterior (sum of log-likelihood and log-prior)
    %   d_thetaR      - The derivative of the log-posterior with respect to
    %                   the random parameters in thetaR
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    if nargout > 1
        % Compute log-likelihood and its derivatives
        [log_likelihood, d_likelihood] = logLikelihood(thetaR, thetaF, varR, varF, C);

        % Compute log-prior and its derivatives
        [log_prior, d_prior] = logPrior(thetaR, thetaP);
        
        % Compute the derivative of the log-posterior
        d_thetaR = d_likelihood + d_prior;
    else
        % Compute log-likelihood without derivatives
        log_likelihood = logLikelihood(thetaR, thetaF, varR, varF, C);

        % Compute log-prior without derivatives
        log_prior = logPrior(thetaR, thetaP);
    end

    % Compute the log-posterior as the sum of log-likelihood and log-prior
    log_posterior = log_likelihood + log_prior;
end