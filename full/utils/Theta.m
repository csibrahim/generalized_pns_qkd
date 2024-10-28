function theta = Theta(phi, thetaP)
    %Theta: This function transforms unconstrained parameters (phi) into 
    %        constrained parameters (theta) using the upper (ub) and lower (lb) bounds 
    %        provided in thetaP. It applies different transformations for bounded and 
    %        unbounded parameters.
    %
    %        For bounded parameters, a sigmoid function is applied to map the values 
    %        to the range [lb, ub]. For unbounded parameters, an exponential function 
    %        is used to ensure positive values, shifted by the lower bound.
    %
    % Inputs:
    %     phi     - Unconstrained parameters to be transformed (size: n x m)
    %     thetaP  - Cell array containing the prior parameters {alphas, betas, ub, lb}
    %
    % Outputs:
    %     theta   - Constrained system parameters, transformed from phi (size: n x m)
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Extract upper and lower bounds from thetaP
    [~,~,ub,lb] = deal(thetaP{:});
    ub = [ub{:}];  % Convert cell array to regular array
    lb = [lb{:}];

    % Identify unbounded (Gamma) and bounded (Beta) parameters
    gs = isinf(ub);  % Unbounded parameters have infinite upper bounds
    bs = ~gs;        % Bounded parameters have finite bounds

    % Compute the width for bounded parameters
    width = ub - lb;

    % Initialize theta with the same size as phi
    theta = zeros(size(phi));

    % Transform unbounded parameters using exponential transformation
    theta(:, gs) = exp(phi(:, gs)) + lb(:, gs);

    % Transform bounded parameters using sigmoid and scaling by width
    theta(:, bs) = sigmoid(phi(:, bs)) .* width(:, bs) + lb(:, bs);
end