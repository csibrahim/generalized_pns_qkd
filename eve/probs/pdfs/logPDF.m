function [log_pdf, d_phiE] = logPDF(C, thetaA, thetaB, phiE, thetaP, factor)
    %logPDF: This function computes the logarithm of the posterior PDF given the 
    %        system parameters and observed counts, incorporating a Jacobian 
    %        adjustment to account for the transformation from the unconstrained 
    %        parameters (phiE) to the constrained parameters (thetaE). It also 
    %        supports optional scaling through the `factor` input and, if requested, 
    %        computes the derivative of the log-PDF with respect to phiE.
    %
    % Inputs:
    %   C        - Array of observed counts for each outcome.
    %   thetaA   - Cell array of system parameters specific to the channel between 
    %              Alice and Bob, containing [lambdas, alpha, dAB].
    %   thetaB   - Cell array of system parameters for Bob's detectors, containing 
    %              [pa0, pa1, pc0, pc1, pd0, pd1, pe].
    %   phiE     - Unconstrained parameters to be transformed into constrained system 
    %              parameters (thetaE).
    %   thetaP   - Cell array containing the prior parameters (alphas, betas, ub, lb).
    %   factor   - (Optional) Scaling factor applied to the log-PDF and its derivative.
    %
    % Outputs:
    %   log_pdf  - The logarithm of the posterior PDF, including the Jacobian adjustment.
    %   d_phiE   - The derivative of the log-PDF with respect to phiE.
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Set scaling factor to 1 if not provided
    if nargin < 6
        factor = 1;
    end

    % Extract prior parameters and compute width for bounded parameters
    [~, ~, ub, lb] = deal(thetaP{:});
    width = ub - lb;

    % Identify bounded (Beta) and unbounded (Gamma) parameters
    gs = isinf(ub);  % Unbounded parameters (Gamma-distributed)
    bs = ~gs;        % Bounded parameters (Beta-distributed)

    % Determine the size of phiE
    [n, m] = size(phiE);

    % Initialize thetaE (transformed system parameters)
    thetaE = zeros(n, m);

    % Apply sigmoid for bounded and exponential for unbounded transformations
    s = sigmoid(phiE(:, bs));
    e = exp(phiE(:, gs));

    % Transform phiE to thetaE
    thetaE(:, bs) = s .* width(bs) + lb(bs);
    thetaE(:, gs) = e + lb(gs);

    % Convert thetaE to cell format for compatibility with other functions
    thetaE = mat2cell(thetaE, n, ones(1, m));

    if nargout > 1
        % Compute log-posterior and its gradient
        [log_post, d_post_d_thetaE] = logPosterior(C, thetaA, thetaB, thetaE, thetaP);

        % Compute log-Jacobian and its gradient
        [log_jacob, d_jacob_d_phiE] = logJacobian(phiE, thetaP);
    
        % Compute the Jacobian of the transformation with respect to phiE
        d_theta_d_phiE = zeros(n, m);    
        d_theta_d_phiE(:, bs) = s .* (1 - s) .* width(bs);  % Sigmoid derivative
        d_theta_d_phiE(:, gs) = e;                          % Exponential derivative

        % Adjust the gradient of the posterior
        d_post_d_phiE = d_post_d_thetaE .* d_theta_d_phiE;
        
        % Combine the gradients for the total derivative of log-PDF
        d_phiE = (d_post_d_phiE + d_jacob_d_phiE) * factor;
    else
        % Compute log-posterior and log-Jacobian without gradients
        log_post = logPosterior(C, thetaA, thetaB, thetaE, thetaP);
        log_jacob = logJacobian(phiE, thetaP);
    end

    % Combine log-posterior and log-Jacobian to compute total log-PDF
    log_pdf = (log_post + log_jacob) * factor;
end