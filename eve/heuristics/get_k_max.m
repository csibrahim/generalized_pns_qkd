function k_max = get_k_max(thetaA, thetaB)
    %get_k_max: Find the maximum photons, k, Eve can intercept such that 
    %           Eve requires a channel transmission probability pEB = 1 to 
    %           remain undetected.
    %
    % Inputs:
    %   thetaA - Alice’s parameter set [lambdas, alpha, dAB]
    %   thetaB - Bob’s parameter set [pa0, pa1, pc0, pc1, pd0, pd1, pe]
    %
    % Outputs:
    %   k_max  - Maximum k requiring pEB = 1 for undetection
    %
    % Copyright (c) 2024 Ibrahim Almosallam
    % Licensed under the MIT License (see LICENSE file for full details).

    % Load optimization options
    opt_options;

    % Extract Alice’s parameters
    [lambdas, alpha, dAB] = thetaA{:};

    % Use the highest intensity in lambdas as the target intensity
    lambda_max = max(lambdas);

    % Redefine thetaA to only include the max intensity
    thetaA = {lambda_max, alpha, dAB};

    % Define the cost function
    f = @(log_k) kCost(log_k, thetaA, thetaB);

    % Use fminunc to find the optimal log(k) that minimizes the
    % KL divergence initialized at log(max_lambda)
    [log_k,~,~,output] = fminunc(f, log(lambda_max), options);
    
    % Convert optimized log(k) to k and add 1 for offset
    k_max = exp(log_k) + 1;
    
    % Warn if optimization fails and fallback to lambda_max or 2
    if output.iterations == 0
        warning("Finding maximum k failed, setting to k_max = max(lambda_max, 2)");
        k_max = max(lambda_max, 2);
    end

    % Ensure k_max is at least 2
    k_max = max(k_max, 2);
    
end

function [KL, d_log_k] = kCost(log_k, thetaA, thetaB)
    %kCost: Computes KL divergence, and its gradient, between:
    %       1) Detector clicks when Eve does not intercept.
    %       2) Detector clicks when Eve intercepts k photons and pEB = 1.
    %
    % Inputs:
    %   log_k   - Log-transformed k for unconstrained optimization.
    %   thetaA  - Alice’s parameters [lambdas, alpha, dAB]
    %   thetaB  - Bob’s parameters [pa0, pa1, pc0, pc1, pd0, pd1, pe]
    %
    % Outputs:
    %   KL      - KL divergence between no-interception and interception cases
    %   d_log_k - Derivative of the KL divergence w.r.t log(k)
    %
    % Copyright (c) 2024 Ibrahim Almosallam
    % Licensed under the MIT License (see LICENSE file for full details).

    % Define Eve's parameters
    dAE = 0;            % Distance between Alice and Eve
    pEB = 1;            % Perfect channel transmission
    k = exp(log_k) + 1; % Convert log(k) back to k
    Delta = 1;          % Eve intercepts all pulses

    % Set Eve's parameter set
    thetaE = {dAE, pEB, k, Delta};

    % Assume Alice and Bob use matching bases
    a = 1; % Alice's basis choice
    b = 1; % Bob's basis choice

    % Compute detection probabilities and gradients
    Ps0 = Pabe(thetaA, thetaB, thetaE, a, b, 0);               % No interception
    [Ps1, d_thetaE1] = Pabe(thetaA, thetaB, thetaE, a, b, 1);  % Intercepting with k photons

    % Compute KL divergence
    KL1 = Ps0 .* log(Ps0); % Expected log likelihood for Ps0
    KL0 = Ps0 .* log(Ps1); % Expected log likelihood for Ps1
    KL = sum(KL1 - KL0);   % KL divergence as sum of differences

    % Gradient of KL divergence w.r.t k
    d_k = d_thetaE1{3};            % Derivative of Ps1 w.r.t k
    d_KL = - (Ps0 ./ Ps1) .* d_k;  % Derivative of KL w.r.t k
    d_log_k = sum(d_KL) * (k - 1); % Derivative of KL w.r.t log(k)

end