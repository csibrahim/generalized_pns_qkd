function [log_beta, d_beta] = logBeta(p, alphas, betas)
    %logBeta: This function computes the logarithm of the Beta distribution 
    %         PDF for given parameters alphas and betas, evaluated at the 
    %         probability p. To avoid numerical issues when p approaches 0 or 
    %         1, the function applies a small scaling to p. If requested, it 
    %         also computes the derivative of the log-Beta function with 
    %         respect to the parameters.
    %
    % Inputs:
    %   p      - Probability value(s) at which the Beta distribution is evaluated
    %   alphas - Shape parameter(s) of the Beta distribution (first parameter)
    %   betas  - Shape parameter(s) of the Beta distribution (second parameter)
    %
    % Outputs:
    %   log_beta - The logarithm of the Beta distribution evaluated at p
    %   d_beta   - The derivative of log_beta with respect to the paramters
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Scale p to avoid numerical issues at boundaries (p = 0 or p = 1)
    p_scaled = (1 - eps) * p + eps / 2;

    % Compute the log terms for alphas and betas
    log_a = (alphas - 1) .* log(p_scaled);
    log_b = (betas  - 1) .* log1p(-p_scaled);
    
    % Compute the logarithm of the Beta distribution (without normalization)
    log_beta = log_a + log_b; % - betaln(alphas, betas);

    % If derivatives are requested
    if nargout > 1
        % Compute the derivative of log-Beta with respect to the parameters
        d_beta = (1 - eps) * ((alphas - 1) ./ p_scaled - (betas - 1) ./ (1 - p_scaled));
    end

end