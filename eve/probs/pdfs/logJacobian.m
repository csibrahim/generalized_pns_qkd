function [log_jacob, grad] = logJacobian(phi, thetaP)
    %logJacobian: This function computes the logarithm of the Jacobian for 
    %             the transformation from the unconstrained variables (phi) 
    %             to the constrained variables (thetaE), with respect to the 
    %             prior parameters. The Jacobian accounts for bounded 
    %             (Beta-distributed) and unbounded (Gamma-distributed) 
    %             parameters. If requested, the function also returns the 
    %             gradient of the log-Jacobian with respect to phi.
    %
    % Inputs:
    %   phi     - Unconstrained parameters to be transformed
    %   thetaP  - Cell array containing the prior parameters (alphas, betas, ub, lb)
    %
    % Outputs:
    %   log_jacob - The logarithm of the Jacobian for the transformation
    %   grad      - The gradient of the log-Jacobian with respect to phi
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Extract upper bounds, lower bounds
    [~, ~, ub, lb] = deal(thetaP{:});
    
    width = ub - lb;  % Compute the width for bounded parameters

    % Identify unbounded (Gamma) and bounded (Beta) parameters
    gs = isinf(ub);  % Unbounded parameters (Gamma-distributed)
    bs = ~gs;        % Bounded parameters (Beta-distributed)

    % Compute sigmoid and its derivative for the bounded parameters
    [s, ds] = sigmoid(phi(:, bs));
    s_scaled = (1 - eps) * s + eps / 2;  % Scale sigmoid to avoid numerical issues

    % Initialize the log-Jacobian matrix (only the diagonal)
    log_jacob = zeros(size(phi));

    % Compute the log-Jacobian for bounded and unbounded parameters
    log_jacob(:, gs) = phi(:,gs);
    log_jacob(:, bs) = log(s_scaled) + log(1 - s_scaled) + log(width(bs));

    % Sum across columns to get the total log-Jacobian
    log_jacob = sum(log_jacob, 2);
    
    % If the gradient is requested
    if nargout > 1
        % Initialize the gradient matrix
        grad = zeros(size(phi));
        
        % Compute the gradient for bounded parameters (Beta-distributed)
        grad(:, bs) = (1 - eps) * ds .* (1 ./ s_scaled - 1 ./ (1 - s_scaled));
        
        % The gradient for unbounded parameters (Gamma-distributed) is 1
        grad(:, gs) = 1;
    end
end