function phi = Phi(theta, thetaP)
    %Phi: This function transforms constrained parameters (theta) into 
    %      unconstrained parameters (phi) using the upper (ub) and lower (lb) bounds 
    %      provided in thetaP. It applies different transformations depending on 
    %      whether the parameters are bounded or unbounded.
    %
    %      For bounded parameters, the inverse sigmoid (logit) function is applied 
    %      to map the values from the constrained space to the unconstrained space. 
    %      For unbounded parameters, the logarithmic function is used to map the values 
    %      from the constrained space.
    %
    % Inputs:
    %     theta   - Constrained system parameters to be transformed (size: n x m)
    %     thetaP  - Cell array containing the prior parameters {alphas, betas, ub, lb}
    %
    % Outputs:
    %     phi     - Unconstrained parameters, transformed from theta (size: n x m)
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Extract upper and lower bounds from thetaP
    [~, ~, ub, lb] = deal(thetaP{:});
    ub = [ub{:}];  % Convert cell array to regular array
    lb = [lb{:}];

    % Identify unbounded (Gamma) and bounded (Beta) parameters
    gs = isinf(ub);  % Unbounded parameters (infinite upper bounds)
    bs = ~gs;        % Bounded parameters (finite bounds)

    % Compute the width for bounded parameters
    width = ub - lb;

    % Initialize phi with the same size as theta
    phi = zeros(size(theta));

    % Transform unbounded parameters using the logarithmic function
    phi(:, gs) = log(theta(:, gs) - lb(gs));  % Logarithmic transformation for unbounded

    % Transform bounded parameters using the inverse sigmoid (logit) function
    phi(:, bs) = logit((theta(:, bs) - lb(bs)) ./ width(bs));  % Logit transformation for bounded
end