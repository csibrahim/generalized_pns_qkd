function [log_likelihood, d_thetaE] = logLikelihood(C, thetaA, thetaB, thetaE)
    %logLikelihood: This function computes the log-likelihood of observing the 
    %               counts C under a multinomial model given the detection 
    %               probabilities P. If requested, it also computes the gradient 
    %               of the log-likelihood with respect to the parameters in thetaE.
    %
    % Inputs:
    %   C        - Array of observed counts for each outcome.
    %   thetaA   - Cell array of system parameters specific to the channel between 
    %              Alice and Bob, containing [lambdas, alpha, dAB].
    %   thetaB   - Cell array of system parameters for Bob's detectors, containing 
    %              [pa0, pa1, pc0, pc1, pd0, pd1, pe].
    %   thetaE   - Cell array of system parameters related to Eve’s eavesdropping, 
    %              containing [dAE, pEB, k, Delta].
    %
    % Outputs:
    %   log_likelihood - The log-likelihood of the observed counts under the 
    %                    multinomial model.
    %   d_thetaE       - Cell array containing the gradient of the log-likelihood 
    %                    with respect to the parameters in thetaE.
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Initialize output cells
    varargout = cell(1, max(nargout, 1));

    % Compute detection probabilities (Ps) and optionally their gradients
    [varargout{:}] = P(thetaA, thetaB, thetaE);
    Ps = varargout{1};  % Extract detection probabilities

    % Add a small epsilon to avoid issues with zero probabilities
    Ps_eps = Ps + eps;

    % Normalize probabilities to ensure they sum to 1 across rows
    S = sum(Ps_eps, 2); 
    Ps_norm = Ps_eps ./ S;

    % Compute the log-likelihood (unnormalized)
    log_likelihood = sum(C .* log(Ps_norm), 2);% -sum(gammaln(C+1),2) + gammaln(sum(C+1,2));

    if (nargout > 1)
        % Compute the gradient of the log-likelihood with respect to thetaE
        d_thetaE = varargout{2};  % Extract gradients of Ps
        
        % Loop over each parameter in thetaE to compute the gradient
        for i = 1:numel(d_thetaE)
            % Derivative of the sum S with respect to the parameter
            dSdx = sum(d_thetaE{i}, 2);

            % Gradient of normalized probabilities
            dPs_norm = (d_thetaE{i} .* S - Ps_eps .* dSdx) ./ S.^2;

            % Compute the gradient of the log-likelihood with respect to the parameter
            d_thetaE{i} = sum(C .* (dPs_norm ./ Ps_norm), 2);
        end

        d_thetaE = [d_thetaE{:}];
    end
end