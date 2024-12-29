function [log_pdf, d_pdf] = logPDF(phiR, thetaF, varR, varF, C, thetaP, factor)
    %logPDF: Computes the logarithm of the posterior probability density 
    %        function (PDF) given system parameters and observed counts. 
    %        Includes a Jacobian adjustment for the transformation from 
    %        unconstrained (phiR) to constrained (thetaR) parameters. 
    %        Optionally calculates the gradient of the log-PDF with respect 
    %        to phiE and supports scaling through an optional factor.
    %
    % Inputs:
    %   phiR    - Unconstrained parameters to be transformed into constrained 
    %             parameters (thetaR) 
    %   thetaF  - Values of parameters listed in varF 
    %   varR    - Cell array of random variable names
    %   varF    - Cell array of fixed variable names
    %   C       - Observed click counts
    %   thetaP  - Prior parameters of random variables [alphas, betas, ub, lb]
    %   factor  - Scaling factor to the log-PDF and its derivative (default: 1)
    %
    % Outputs:
    %   log_pdf  - The logarithm of the posterior PDF, including the Jacobian adjustment
    %   d_pdf    - The Derivative of the log-PDF with respect to phiR
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Set the scaling factor to 1 if not provided
    if nargin < 7
        factor = 1;
    end

    % Extract prior parameters: upper bounds, lower bounds, and their differences
    [~, ~, ub, lb] = deal(thetaP{:});

    % Flatten parameters
    ub = [ub{:}];
    lb = [lb{:}];

    % Compute width for bounded parameters
    width = ub - lb;

    % Identify semi-bounded (Gamma) and bounded (Beta) parameters
    gs = isinf(ub);  % Semi-bounded parameters (Gamma-distributed)
    bs = ~gs;        % Bounded parameters (Beta-distributed)

    % Determine the size of phiR
    [n, m] = size(phiR);

    % Initialize thetaR
    thetaR = zeros(n, m);

    % Apply transformations: sigmoid for bounded, exponential for semi-bounded
    s = sigmoid(phiR(:, bs)); % Sigmoid transformation for bounded parameters
    e = exp(phiR(:, gs));     % Exponential transformation for semi-bounded parameters

    % Transform phiR to thetaR
    thetaR(:, bs) = s .* width(bs) + lb(bs);
    thetaR(:, gs) = e + lb(gs);

    % Convert thetaR from matrix to cell format
    thetaR = theta2cell(thetaR, varR);

    % If derivatives are requested
    if nargout > 1

        % Compute log-posterior and its derivatives
        [log_post, d_post_d_thetaR] = logPosterior(thetaR, thetaF, varR, varF, C, thetaP);

        % Compute log-Jacobian and its derivatives
        [log_jacob, d_jacob_d_phiR] = logJacobian(phiR, thetaP);
        
        % Compute the derivative of thetaR w.r.t. phiR
        d_thetaR_d_phiR = zeros(n, m);    
        d_thetaR_d_phiR(:, bs) = s .* (1 - s) .* width(bs);  % Sigmoid derivative for Beta-distributed parameters
        d_thetaR_d_phiR(:, gs) = e;                          % Exponential derivative for Gamma-distributed parameters

        % Apply the chain rule
        d_post_d_thetaR = d_post_d_thetaR .* d_thetaR_d_phiR;

        % Compute the total derivative of the log-PDF
        d_pdf = (d_post_d_thetaR + d_jacob_d_phiR) * factor;
    else

        % Compute log-posterior without derivatives
        log_post = logPosterior(thetaR, thetaF, varR, varF, C, thetaP);

        % Compute log-Jacobian without derivative
        log_jacob = logJacobian(phiR, thetaP);
    end

    % Compute the final log-PDF
    log_pdf = (log_post + log_jacob) * factor;

end