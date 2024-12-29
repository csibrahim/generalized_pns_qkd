function [log_beta, d_beta] = logBeta(p, alphas, betas)
    %logBeta: Computes the logarithm of the Beta distribution for a given 
    %         probability p and parameters alphas and betas. Optionally
    %         calculates the derivative of log-Beta w.r.t the parameters.
    %
    % Inputs:
    %   p      - Probability value(s) where the Beta distribution is evaluated
    %   alphas - Shape parameter(s) of the Beta distribution
    %   betas  - Shape parameter(s) of the Beta distribution
    %
    % Outputs:
    %   log_beta - Logarithm of the Beta distribution evaluated at p
    %   d_beta   - Derivative of log_beta w.r.t. the parameters:
    %              - d_beta(:, 1): Derivative w.r.t. alphas
    %              - d_beta(:, 2): Derivative w.r.t. betas
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Scale p to ensure numerical stability near boundaries (p = 0 or p = 1)
    p_scaled = (1 - eps) * p + eps / 2;

    % Compute the logarithm of the Beta distribution (without normalization)
    log_beta = (alphas - 1) .* log(p_scaled) + (betas - 1) .* log1p(-p_scaled); % - betaln(alpha, beta)

    % Compute derivatives if requested
    if nargout > 1
        % Derivative of log_beta w.r.t. alphas and betas
        d_beta = (1 - eps) * ((alphas - 1) ./ p_scaled - (betas - 1) ./ (1 - p_scaled));
    end
end
