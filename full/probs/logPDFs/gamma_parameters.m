function [alpha, beta, ub, lb] = gamma_parameters(E, V)
    %gamma_parameters: Computes the alpha and beta parameters for a 
    %                  Gamma distribution given the expected value (E) and 
    %                  the variance (V).
    %
    % Inputs:
    %     E - The expected value of the Gamma distribution
    %     V - The variance of the Gamma distribution
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
    alpha = E.^2 ./ V;

    % Compute the beta (rate) parameter
    beta = E ./ V;

    % Set the upper and lower bounds for the Gamma distribution
    ub = inf(size(E));    % Upper bound: infinity
    lb = zeros(size(E));  % Lower bound: 0

end