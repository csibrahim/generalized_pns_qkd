function [log_pdf, d_phiE] = logPDF(C, thetaA, thetaB, phiE, thetaP, factor)
    %logPDF: Computes the logarithm of the posterior probability density 
    %        function (PDF) given system parameters and observed counts. 
    %        Includes a Jacobian adjustment for the transformation from 
    %        unconstrained (phiE) to constrained (thetaE) parameters. 
    %        Optionally calculates the gradient of the log-PDF with respect 
    %        to phiE and supports scaling through an optional factor.
    %
    % Inputs:
    %   C        - Array of observed counts for each outcome
    %   thetaA   - Alice’s parameter set [lambdas, alpha, dAB]
    %   thetaB   - Bob’s parameter set [pa0, pa1, pc0, pc1, pd0, pd1, pe]
    %   phiE     - Unconstrained parameters to be transformed into constrained 
    %              parameters (thetaE)
    %   thetaP   - Priors for Eve’s parameters [alphas, betas, ub, lb]
    %   factor   - Scaling factor to the log-PDF and its derivative (default: 1)
    %
    % Outputs:
    %   log_pdf  - The logarithm of the posterior PDF, including the Jacobian adjustment
    %   d_pdf    - The Derivative of the log-PDF with respect to phiE
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Set the scaling factor to 1 if not provided
    if nargin < 6
        factor = 1;
    end

    % Extract prior parameters
    [~, ~, ub, lb] = deal(thetaP{:});

    % Compute width for bounded parameters
    width = ub - lb;

    % Identify semi-bounded (Gamma) and bounded (Beta) parameters
    gs = isinf(ub);  % Semi-bounded parameters (Gamma-distributed)
    bs = ~gs;        % Bounded parameters (Beta-distributed)

    % Determine the size of phiE
    [n, m] = size(phiE);

    % Initialize thetaE
    thetaE = zeros(n, m);

    % Apply transformations: sigmoid for bounded, exponential for semi-bounded
    s = sigmoid(phiE(:, bs));  % Sigmoid transformation for bounded parameters
    e = exp(phiE(:, gs));      % Exponential transformation for semi-bounded parameters

    % Transform phiE to thetaE
    thetaE(:, bs) = s .* width(bs) + lb(bs);
    thetaE(:, gs) = e + lb(gs);

    % Convert thetaE from matrix to cell format
    thetaE = mat2cell(thetaE, n, ones(1, m));

    % If derivatives are requested
    if nargout > 1

        % Compute log-posterior and its derivatives
        [log_post, d_post_d_thetaE] = logPosterior(C, thetaA, thetaB, thetaE, thetaP);

        % Compute log-Jacobian and its derivatives
        [log_jacob, d_jacob_d_phiE] = logJacobian(phiE, thetaP);

        % Compute the derivative of thetaE w.r.t. phiE
        d_thetaE_d_phiE = zeros(n, m);
        d_thetaE_d_phiE(:, bs) = s .* (1 - s) .* width(bs);  % Sigmoid derivative for Beta-distributed parameters
        d_thetaE_d_phiE(:, gs) = e;                          % Exponential derivative for Gamma-distributed parameters

        % Apply the chain rule
        d_post_d_phiE = d_post_d_thetaE .* d_thetaE_d_phiE;

        % Combine gradients to compute total derivative of the log-PDF
        d_phiE = (d_post_d_phiE + d_jacob_d_phiE) * factor;
    else
        % Compute log-posterior without derivatives
        log_post = logPosterior(C, thetaA, thetaB, thetaE, thetaP);

        % Compute log-Jacobian without derivative
        log_jacob = logJacobian(phiE, thetaP);
    end

    % Compute the final log-PDF
    log_pdf = (log_post + log_jacob) * factor;
end
