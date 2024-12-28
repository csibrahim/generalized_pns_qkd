function [log_likelihood, d_thetaR] = logLikelihood(thetaR, thetaF, varR, varF, C)
    %logLikelihood: Computes the log-likelihood of observed counts C under a 
    %               multinomial distribution given detection probabilities.
    %               Optionally calculates the gradient of the log-likelihood
    %               w.r.t. the parameters in thetaR.
    %
    % Inputs:
    %   thetaR  - Values of parameters listed in varR
    %   thetaF  - Values of parameters listed in varF
    %   varR    - Cell array of random variable names
    %   varF    - Cell array of fixed variable names
    %   C       - Observed click counts
    %
    % Outputs:
    %   log_likelihood - Log-likelihood of observed counts under the
    %                    multinomial distribution
    %   d_thetaR       - Derivatives of the log-likelihood w.r.t. parameters 
    %                    in thetaR
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Initialize output cells
    varargout = cell(1, max(nargout, 1));

    % Compute detection probabilities (Ps) and optionally their gradients
    [varargout{:}] = P(thetaR, thetaF, varR, varF);
    Ps = varargout{1};  % Extract detection probabilities

    % Add a small epsilon to avoid issues with zero probabilities
    Ps_eps = Ps + eps;

    % Normalize probabilities to ensure they sum to 1 across rows
    S = sum(Ps_eps, 2); 
    Ps_norm = Ps_eps ./ S;

    % Compute the log-likelihood based (unnormalized)
    log_likelihood = sum(C .* log(Ps_norm), 2); % - sum(gammaln(C+1),2) + gammaln(sum(C+1,2));

    if nargout > 1
        % Compute the derivative of the log-likelihood w.r.t. thetaR
        d_thetaR = varargout{2};  % Extract derivatives of Ps
        
        % Compute the derivative w.r.t. each parameter in thetaR
        for i = 1:numel(d_thetaR)
            % Derivative of the sum S w.r.t. the parameter
            dSdx = sum(d_thetaR{i}, 2);

            % Gradient of normalized probabilities
            dPs_norm = (d_thetaR{i} .* S - Ps_eps .* dSdx) ./ S.^2;

            % Derivative of the log-likelihood w.r.t. thetaR{i}
            d_thetaR{i} = sum(C .* (dPs_norm ./ Ps_norm), 2);
        end

        % Combine derivatives into a single array
        d_thetaR = [d_thetaR{:}];
    end
end