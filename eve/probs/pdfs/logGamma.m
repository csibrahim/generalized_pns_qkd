function [log_gamm, d_gamma] = logGamma(x, alpha, beta)
    %logGamma: This function computes the logarithm of the Gamma distribution 
    %          PDF for a given input x and parameters alpha and beta, with a 
    %          small adjustment to x to avoid numerical issues near zero. If 
    %          requested, the function also computes the derivative of the 
    %          log-Gamma with respect to the parameters.
    %
    % Inputs:
    %   x      - The input value(s) at which the Gamma distribution is evaluated
    %   alpha  - Shape parameter of the Gamma distribution
    %   beta   - Rate parameter of the Gamma distribution
    %
    % Outputs:
    %   log_gamm - The logarithm of the Gamma distribution evaluated at x
    %   d_gamma  - The derivative of log_gamm with respect to the parameters
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Adjust x to avoid numerical issues (ensures x > 0)
    x = x + eps;
    
    % Compute the logarithm of the Gamma distribution (without normalization)
    log_gamm = (alpha - 1) .* log(x + (alpha == 1)) - beta .* x; % + alpha .* log(beta) - gammaln(alpha);

    % If derivatives are requested
    if nargout > 1
        % Compute the derivative of log-Gamma with respect to the
        % parameters.
        d_gamma = (alpha - 1) ./ (x + (alpha == 1)) - beta;
    end

end
