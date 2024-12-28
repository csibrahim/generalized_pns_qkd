function [alpha, beta, ub, lb] = beta_parameters(E, V)
    %beta_parameters: Computes the alpha and beta parameters of a Beta 
    %                 distribution based on the given expected value (E) 
    %                 and variance (V).
    %
    % Inputs:
    %   E - Expected value of the Beta distribution (0 ≤ E ≤ 1)
    %   V - Variance of the Beta distribution (0 ≤ V ≤ E * (1 - E))
    %
    % Outputs:
    %   alpha  - Alpha parameter of the Beta distribution
    %   beta   - Beta parameter of the Beta distribution
    %   ub     - Upper bounds for the Beta distribution (always 1)
    %   lb     - Lower bounds for the Beta distribution (always 0)
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Validate the variance condition
    if any(V > E .* (1 - E))
        error('Invalid variance: V exceeds theoretical limit (V > E * (1 - E)) for the Beta distribution');
    end

    % Compute the concentration parameter
    concentration = E .* (1 - E) ./ V - 1;

    % Compute alpha and beta parameters from the concentration
    alpha = E .* concentration;
    beta = (1 - E) .* concentration;

    % Set the bounds for the Beta distribution (always [0, 1])
    ub = ones(size(E));  % Upper bound: 1
    lb = zeros(size(E)); % Lower bound: 0
end
