function theta = Theta(phi, thetaP)
    %Theta: Transforms unconstrained parameters (phi) into constrained
    %       parameters (theta) using upper (ub) and lower (lb) bounds from
    %       thetaP. Uses sigmoid for bounded parameters and exponential for
    %       semi-bounded parameters.
    %
    % Inputs:
    %     phi     - Unconstrained parameters (size: n x m)
    %     thetaP  - Cell array with prior parameters {alphas, betas, ub, lb}
    %
    % Outputs:
    %     theta   - Constrained parameters transformed from phi (size: n x m)
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Extract upper and lower bounds from thetaP
    [~,~,ub,lb] = deal(thetaP{:});
    
    % Identify parameter types
    gs = isinf(ub);  % Semi-bounded parameters have infinite upper bounds
    bs = ~gs;        % Bounded parameters have finite bounds

    % Calculate width for bounded parameters
    width = ub - lb;

    % Initialize theta with the same size as phi
    theta = zeros(size(phi));

    % Apply exponential transformation for semi-bounded parameters
    theta(:, gs) = exp(phi(:, gs)) + lb(:, gs);

    % Apply sigmoid transformation for bounded parameters
    theta(:, bs) = sigmoid(phi(:, bs)) .* width(:, bs) + lb(:, bs);
end
