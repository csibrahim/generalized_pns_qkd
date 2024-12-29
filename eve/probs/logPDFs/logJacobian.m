function [log_jacob, grad] = logJacobian(phiE, thetaP)
    %logJacobian: Computes the logarithm of the Jacobian for transforming 
    %             unconstrained variables (phiE) to constrained variables 
    %             (thetaE) based on the prior parameters. The Jacobian 
    %             accounts for bounded (Beta-distributed) and semi-bounded 
    %             (Gamma-distributed) parameters. Optionally calculates the
    %             derivative of the log-Jacobian w.r.t. phi.
    %
    % Inputs:
    %   phiE    - Unconstrained parameters to be transformed
    %   thetaP  - Priors for Eve’s parameters [alphas, betas, ub, lb]
    %
    % Outputs:
    %   log_jacob - Logarithm of the Jacobian determinant
    %   grad      - Derivative of the log-Jacobian with respect to phiE
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Extract prior bounds
    [~, ~, ub, lb] = deal(thetaP{:});
    
    % Compute the width for bounded parameters
    width = ub - lb;

    % Identify parameter types
    gs = isinf(ub);  % Semi-bounded (Gamma-distributed)
    bs = ~gs;        % Bounded (Beta-distributed)

    % Compute sigmoid and its derivative for bounded parameters
    [s, ds] = sigmoid(phiE(:, bs));
    s_scaled = (1 - eps) * s + eps / 2;  % Scale sigmoid to avoid numerical issues

    % Initialize the log-Jacobian
    log_jacob = zeros(size(phiE));

    % Compute the log-Jacobian for bounded and semi-bounded parameters
    log_jacob(:, gs) = phiE(:, gs); % Semi-bounded parameters
    log_jacob(:, bs) = log(s_scaled) + log(1 - s_scaled) + log(width(bs)); % Bounded parameters

    % Sum across columns to compute the log-determinant
    log_jacob = sum(log_jacob, 2);

    % If derivatives are requested
    if nargout > 1
        % Initialize the derivative matrix
        grad = zeros(size(phiE));

        % Derivatives for bounded parameters
        grad(:, bs) = (1 - eps) * ds .* (1 ./ s_scaled - 1 ./ (1 - s_scaled));

        % Derivatives for semi-bounded parameters
        grad(:, gs) = 1; 
    end
end
