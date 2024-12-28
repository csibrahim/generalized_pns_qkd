function varargout = PNS(lambdas, pc, pd, k, outs)
    %PNS: Computes the detection probability and their derivatives with 
    %     respect to the input parameters in a Photon-Number Splitting (PNS) 
    %     attack. This assumes Alice, Bob, and Eve are directly connected, 
    %     with Eve intercepting a fixed number of photons (k).
    %
    % Inputs:
    %   lambdas - Array of average photon numbers (mean intensity per pulse)
    %   pc      - Detection efficiency of Bob's detector (0 ≤ pc ≤ 1)
    %   pd      - Dark count probability of Bob's detector (0 ≤ pd ≤ 1)
    %   k       - Number of photons intercepted by Eve
    %   outs    - Logical array (1x5) specifying which outputs to compute:
    %             - outs(1): Detection probability (Ps)
    %             - outs(2): Partial derivative w.r.t. lambda (d_lambda)
    %             - outs(3): Partial derivative w.r.t. pc (d_pc)
    %             - outs(4): Partial derivative w.r.t. pd (d_pd)
    %             - outs(5): Partial derivative w.r.t. k (d_k)
    %
    % Outputs:
    %   Ps        - Detection probabilities
    %   d_lambda  - Partial derivative w.r.t. lambda (mean intensity)
    %   d_pc      - Partial derivative w.r.t. pc (detection efficiency)
    %   d_pd      - Partial derivative w.r.t. pd (dark count probability)
    %   d_k       - Partial derivative w.r.t. k (photons intercepted by Eve)
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Initialize outputs as a cell array
    varargout = cell(1, nargout);

    % Check for optional input `outs` and ensure proper length and consistency
    if nargin < 5
        outs = false(1, 5); % Default: only Ps is computed
        outs(1:nargout) = true;
    else
        if length(outs) ~= 5
            error("`outs` should be a logical array of length 5.");
        end
        if sum(outs) ~= nargout
            error("Mismatch between the number of requested and assigned outputs.");
        end
    end

    % Adjust dimensions of inputs if necessary
    [nl, ml] = size(lambdas);
    [np, mp] = size(pc);

    % Replicate inputs if needed to handle vectorized operations
    if nl == 1 && np > 1
        lambdas = repmat(lambdas, np, 1);
    end
    if mp == 1 && ml > 1
        pc = repmat(pc, 1, ml);
        pd = repmat(pd, 1, ml);
        k = repmat(k, 1, ml);
    end
    
    % Precompute intermediate values
    qd = 1 - pd; % Complement of dark count probability
    qc = 1 - pc; % Complement of detection efficiency
    
    log_qd = log(qd); 
    log_qc = log(qc); 
    log_pc = log(pc);

    qc_lambdas = qc .* lambdas;
    pc_lambdas = pc .* lambdas;
    
    % Precompute terms involving the incomplete gamma function
    log_gammainc_qc_lambdas_k = log(gammainc(qc_lambdas, k));
    log_qd_k_log_qc_pc_lambdas = log_qd - k .* log(qc + (k == 0)) - pc_lambdas;
    e1 = log_qd_k_log_qc_pc_lambdas + log_gammainc_qc_lambdas_k; % First term

    % Precompute the second term if needed
    if outs(1) || outs(4)
        e2 = log_qd + log(gammainc(lambdas, k, 'upper'));
    end

    % Precompute derivative terms if needed
    if outs(2) || outs(3)
        ln_dgdx_lambdas_k = ln_dgdx(lambdas, k);
    end

    % Compute the detection probabilities
    if outs(1)
        Ps = 1 - exp(e1) - exp(e2); % Combine terms to compute detection probabilities
        varargout{1} = Ps;
    end

    % Compute the partial derivative w.r.t. lambda
    if outs(2)
        dL1 = exp(e1 + log_pc) - exp(log_qd + ln_dgdx_lambdas_k);
        dL2 = exp(log_qd + ln_dgdx_lambdas_k);
        d_lambda = dL1 + dL2;
        varargout{sum(outs(1:2))} = d_lambda;
    end

    % Compute the partial derivative w.r.t. pc
    if outs(3)
        log_lambdas = log(lambdas);
        d_pc = -exp(e1 + log(k) - log_qc) + exp(e1 + log_lambdas) + ...
               exp(log_qd - log_qc + ln_dgdx_lambdas_k + log_lambdas);
        varargout{sum(outs(1:3))} = d_pc;
    end

    % Compute the partial derivative w.r.t. pd
    if outs(4)
        d_pd = exp(e1 - log_qd) + exp(e2 - log_qd);
        varargout{sum(outs(1:4))} = d_pd;
    end

    % Compute the partial derivative w.r.t. k
    if outs(5)
        dgdk_qc_lambdas_k = dgdk_numeric(qc_lambdas, k);
        dgdk_lambdas_k = dgdk_numeric(lambdas, k);

        dk1 = exp(e1 + log(abs(log_qc))) .* sign(log_qc) - ...
              exp(log_qd_k_log_qc_pc_lambdas + log(abs(dgdk_qc_lambdas_k))) .* sign(dgdk_qc_lambdas_k);
        dk2 = exp(log_qd + log(abs(dgdk_lambdas_k))) .* sign(dgdk_lambdas_k);

        d_k = dk1 + dk2;
        varargout{sum(outs)} = d_k;
    end
end


% Function to compute the derivative of the regularized lower incomplete gamma function w.r.t. x
function g = ln_dgdx(x, k)
    %ln_dgdx: Computes the derivative of the regularized lower incomplete 
    %         gamma function with respect to x. 

    g = (k - 1) .* log(x + (k == 1)) - x - gammaln(k); % Main expression
    g(k == 0) = -inf;       % Handle special case when k = 0
    g(k == 0 & x == 0) = 0; % Special case handling when x = 0
end

% Function to compute the numerical derivative of the regularized lower incomplete gamma function w.r.t. k
function g = dgdk_numeric(x, k)
    %dgdk_numeric: Computes the numerical derivative of the regularized lower
    %              incomplete gamma function with respect to k using finite 
    %              differences.

    delta = sqrt(eps) * max(k, 1);
    plus_k = k + delta;
    minus_k = max(k - delta, 0);
    g = (gammainc(x, plus_k) - gammainc(x, minus_k)) ./ (plus_k - minus_k);
    g(k == 0 & x == 0) = -realmax;
end

% Function to compute the analytical derivative of the regularized lower incomplete gamma function w.r.t. k
function g = dgdk_analytical(x, k)
    %dgdk_analytical: Computes the analytical derivative of the regularized 
    %                 lower incomplete gamma function with respect to k 
    %                 using integral evaluation. Returns 0 for x = 0.

    integrand = @(t) t.^(k - 1) .* exp(-t) .* log(t);
    if x == 0
        g = 0;
    else
        g = (1 / gamma(k)) * integral(integrand, 0, x) - gammainc(x, k) * psi(k);
    end
end