function [log_pdf, d_pdf] = logPDF(phiR, thetaF, varR, varF, C, thetaP, factor)
    %logPDF: This function computes the logarithm of the posterior PDF given the 
    %        system parameters and observed counts, incorporating a Jacobian 
    %        adjustment to account for the transformation from the unconstrained 
    %        variables (phiR) to the constrained variables (thetaR). The function 
    %        also supports optional scaling via the `factor` input and computes 
    %        the derivative of the log-PDF with respect to phiR if requested.
    %
    % Inputs:
    %   phiR     - Unconstrained parameters to be transformed into system parameters
    %   thetaF   - Array of fixed system parameters
    %   varR     - Cell array of variable names corresponding to the random parameters
    %   varF     - Cell array of variable names corresponding to the fixed parameters
    %   C        - Array of observed counts for each outcome
    %   thetaP   - Cell array containing the prior parameters (alphas, betas, ub, lb)
    %   factor   - (Optional) Scaling factor applied to the log-PDF and its derivative
    %
    % Outputs:
    %   log_pdf  - The logarithm of the posterior PDF, including the Jacobian adjustment
    %   d_pdf    - The derivative of the log-PDF with respect to phiR
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Set the scaling factor to 1 if not provided
    if nargin < 7
        factor = 1;
    end

    % Extract prior parameters: upper bounds, lower bounds, and their differences
    [~, ~, ub, lb] = deal(thetaP{:});
    ub = [ub{:}];
    lb = [lb{:}];
    width = ub - lb;

    % Identify unbounded (Gamma) and bounded (Beta) parameters
    gs = isinf(ub);  % Unbounded parameters (Gamma-distributed)
    bs = ~gs;        % Bounded parameters (Beta-distributed)

    % Determine the size of the unconstrained input phiR
    [n, m] = size(phiR);

    % Initialize thetaR (transformed system parameters)
    thetaR = zeros(n, m);

    % Apply sigmoid to the bounded parameters and exponential to the unbounded ones
    s = sigmoid(phiR(:, bs));
    e = exp(phiR(:, gs));

    % Transform phiR to thetaR
    thetaR(:, bs) = s .* width(bs) + lb(bs);
    thetaR(:, gs) = e + lb(gs);

    % Convert thetaR from matrix to cell format
    thetaR = theta2cell(thetaR, varR);

    % Initialize output variables for posterior and Jacobian
    varargout_post = cell(1, max(nargout, 1));
    varargout_jacob = cell(1, max(nargout, 1));

    % Compute the log-posterior and the Jacobian
    [varargout_post{:}] = logPosterior(thetaR, thetaF, varR, varF, C, thetaP);
    [varargout_jacob{:}] = logJacobian(phiR, thetaP);

    % Extract log-posterior and log-Jacobian values
    log_post = varargout_post{1};
    log_jacob = varargout_jacob{1};

    % Compute the total log-PDF
    log_pdf = log_post + log_jacob;
    log_pdf = factor * log_pdf;

    % If derivatives are requested
    if nargout > 1
        % Compute the derivative of the posterior and the Jacobian
        d_post = varargout_post{2};
        d_jacob = varargout_jacob{1};

        % Compute the derivative of phiR with respect to the parameters
        d_phiR = zeros(n, m);    
        d_phiR(:, bs) = s .* (1 - s) .* width(bs);  % Sigmoid derivative for Beta-distributed parameters
        d_phiR(:, gs) = e;                          % Exponential derivative for Gamma-distributed parameters

        % Adjust the derivative of the posterior
        d_post = d_post .* d_phiR;

        % Compute the total derivative of the log-PDF
        d_pdf = d_post + d_jacob;
        d_pdf = factor * d_pdf;
    end
end