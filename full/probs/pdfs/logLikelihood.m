function [log_likelihood, d_thetaR] = logLikelihood(thetaR, thetaF, varR, varF, C)
    %logLikelihood: This function computes the log-likelihood of observing the 
    %               counts C under a multinomial model given the detection 
    %               probabilities P. If requested, it also computes the gradient 
    %               of the log-likelihood with respect to the random parameters 
    %               in thetaR.
    %
    % Inputs:
    %   thetaR  - Array of system parameters treated as random variables
    %   thetaF  - Array of system parameters treated as fixed variables
    %   varR    - Cell array of variable names corresponding to the random 
    %                  parameters
    %   varF    - Cell array of variable names corresponding to the fixed 
    %                  parameters
    %   C       - Array of observed counts for each outcome
    %
    % Outputs:
    %   log_likelihood - The log-likelihood of the observed counts under the 
    %                    multinomial model
    %   d_thetaR       - (Optional) Gradient of the log-likelihood with respect 
    %                    to the random variables in thetaR
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    varargout = cell(1, max(nargout, 1));

    % Compute detection probabilities (Ps) and optionally their gradients
    [varargout{:}] = P(thetaR, thetaF, varR, varF);
    Ps = varargout{1};  % Extract detection probabilities

    % To avoid issues with zero probabilities, we add a small epsilon value
    Ps_eps = Ps + eps;

    % Normalize the probabilities to ensure they sum to 1 across rows
    S = sum(Ps_eps, 2); 
    Ps_norm = Ps_eps ./ S;

    % Compute the log-likelihood based on the normalized probabilities
    log_likelihood = sum(C .* log(Ps_norm), 2);

    if nargout > 1
        % Compute the gradient of the log-likelihood with respect to thetaR
        d_thetaR = varargout{2};  % Extract gradients of Ps
        
        % Loop over each random parameter to compute the gradient of log-likelihood
        for i = 1:numel(d_thetaR)
            % Derivative of the sum S w.r.t. the parameter
            dSdx = sum(d_thetaR{i}, 2);

            % Gradient of normalized probabilities
            dPs_norm = (d_thetaR{i} .* S - Ps_eps .* dSdx) ./ S.^2;

            % Compute the gradient of the log-likelihood w.r.t. the parameter
            d_thetaR{i} = sum(C .* (dPs_norm ./ Ps_norm), 2);
        end

        % Combine the gradients into a single output matrix
        d_thetaR = [d_thetaR{:}];
    end
end