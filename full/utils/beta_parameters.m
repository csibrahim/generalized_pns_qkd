function [alpha, beta, ub, lb] = beta_parameters(mus, sigmas)
    %beta_parameters: This function computes the alpha and beta parameters for the 
    %                  Beta distribution, given the means (mus) and standard deviations 
    %                  (sigmas). It also returns the upper and lower bounds for the parameters.
    %
    % Inputs:
    %     mus    - Mean values (mu) for the Beta distribution
    %     sigmas - Standard deviations (sigma) corresponding to each mean
    %
    % Outputs:
    %     alpha  - Alpha parameter for the Beta distribution
    %     beta   - Beta parameter for the Beta distribution
    %     ub     - Upper bounds (always 1 for Beta distribution)
    %     lb     - Lower bounds (always 0 for Beta distribution)
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Ensure that the variance condition is satisfied
    if any(sigmas.^2 > mus .* (1 - mus))
        error('Invalid variance: Variance V exceeds theoretical limit (V > mu * (1 - mu)) for the Beta distribution');
    end

    % Compute the concentration parameter
    concentration = mus .* (1 - mus) ./ sigmas.^2 - 1;

    % Compute alpha and beta parameters from the concentration
    alpha = mus .* concentration;
    beta = (1 - mus) .* concentration;

    % Set the upper and lower bounds for the Beta distribution (always [0, 1])
    ub = ones(size(mus));  % Upper bound: 1
    lb = zeros(size(mus)); % Lower bound: 0

end