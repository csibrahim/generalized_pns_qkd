function [alpha, beta, ub, lb] = gamma_parameters(mus, sigmas)
    %gamma_parameters: Computes the alpha and beta parameters for the 
    %                   Gamma distribution, given the means (mus) and standard deviations 
    %                   (sigmas). It also returns the upper and lower bounds for the parameters.
    %
    % Inputs:
    %     mus    - Mean values (mu) for the Gamma distribution
    %     sigmas - Standard deviations (sigma) corresponding to each mean
    %
    % Outputs:
    %     alpha  - Alpha parameter for the Gamma distribution (shape)
    %     beta   - Beta parameter for the Gamma distribution (rate)
    %     ub     - Upper bounds (always infinity for the Gamma distribution)
    %     lb     - Lower bounds (always 0 for the Gamma distribution)
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Compute the alpha (shape) parameter
    alpha = mus.^2 ./ sigmas.^2;

    % Compute the beta (rate) parameter
    beta = mus ./ sigmas.^2;

    % Set the upper and lower bounds for the Gamma distribution
    ub = inf(size(mus));    % Upper bound: infinity
    lb = zeros(size(mus));  % Lower bound: 0

end