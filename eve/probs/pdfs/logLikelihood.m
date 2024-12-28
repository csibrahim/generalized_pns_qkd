function [log_likelihood, d_thetaE] = logLikelihood(C, thetaA, thetaB, thetaE)
    %logLikelihood: Computes the log-likelihood of observed counts C under a 
    %               multinomial distribution given detection probabilities.
    %               Optionally calculates the gradient of the log-likelihood
    %               w.r.t. the parameters in thetaE.
    %
    % Inputs:
    %   C       - Observed click counts
    %   thetaA  - Alice’s parameter set [lambdas, alpha, dAB]
    %   thetaB  - Bob’s parameter set [pa0, pa1, pc0, pc1, pd0, pd1, pe]
    %   thetaE  - Eve’s parameter set [dAE, pEB, k, Delta]
    %
    % Outputs:
    %   log_likelihood - Log-likelihood of observed counts under the
    %                    multinomial distribution
    %   d_thetaE       - Derivatives of the log-likelihood w.r.t. parameters 
    %                    in thetaE
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

    % Compute the log-likelihood based (unnormalized)
    log_likelihood = sum(C .* log(Ps_norm), 2); % - sum(gammaln(C+1),2) + gammaln(sum(C+1,2));

    if nargout > 1
        % Compute the derivative of the log-likelihood w.r.t. thetaE
        d_thetaE = varargout{2};  % Extract derivatives of Ps
        
        % Compute the derivative w.r.t. each parameter in thetaE
        for i = 1:numel(d_thetaE)
            % Derivative of the sum S w.r.t. the parameter
            dSdx = sum(d_thetaE{i}, 2);

            % Gradient of normalized probabilities
            dPs_norm = (d_thetaE{i} .* S - Ps_eps .* dSdx) ./ S.^2;

            % Derivative of the log-likelihood w.r.t. thetaE{i}
            d_thetaE{i} = sum(C .* (dPs_norm ./ Ps_norm), 2);
        end

        % Combine derivatives into a single array
        d_thetaE = [d_thetaE{:}];
    end
end
