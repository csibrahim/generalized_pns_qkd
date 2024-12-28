function phi = Phi(theta, thetaP)
    %Phi: Transforms constrained parameters (theta) into unconstrained
    %     parameters (phi) using upper (ub) and lower (lb) bounds from
    %     thetaP. Applies logit transformation for bounded parameters and
    %     logarithmic transformation for semi-bounded parameters.
    %
    % Inputs:
    %     theta   - Constrained system parameters (size: n x m)
    %     thetaP  - Cell array with prior parameters {alphas, betas, ub, lb}
    %
    % Outputs:
    %     phi     - Unconstrained parameters transformed from theta (size: n x m)
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Extract bounds from thetaP
    [~, ~, ub, lb] = deal(thetaP{:});

    % Convert cell arrays to regular arrays
    ub = [ub{:}];
    lb = [lb{:}];

    % Identify parameter types
    gs = isinf(ub);  % Semi-bounded parameters (infinite upper bounds)
    bs = ~gs;        % Bounded parameters (finite bounds)

    % Calculate width for bounded parameters
    width = ub - lb;

    % Initialize phi with the same size as theta
    phi = zeros(size(theta));

    % Apply logarithmic transformation for semi-bounded parameters
    phi(:, gs) = log(theta(:, gs) - lb(:, gs));

    % Apply logit transformation for bounded parameters
    phi(:, bs) = logit((theta(:, bs) - lb(:, bs)) ./ width(:, bs));
end
