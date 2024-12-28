function pEB = get_pEB(thetaA, thetaB, dAE, k)
    % get_pEB - Optimizes Eve's transmission probability (pEB) to minimize 
    %           the KL divergence between Bob's detection distributions for:  
    %           1. Eve intercepts k photons and transmits the rest with pEB  
    %           2. Direct transmission without interception
    %
    % Inputs:
    %   thetaA - Alice’s parameters [lambdas, alpha, dAB]
    %   thetaB - Bob’s parameters [pa0, pa1, pc0, pc1, pd0, pd1, pe]
    %   dAE    - Distance from Alice to Eve
    %   k      - Number of photons intercepted by Eve
    %
    % Output:
    %   pEB    - Optimal transmission probability minimizing KL divergence
    %
    % Copyright (c) 2024 Ibrahim Almosallam
    % Licensed under the MIT License (see LICENSE file for full details).
    
    % Set optimization options
    opt_options;
    
    % Define the cost function
    f = @(logit_pEB) pEBcost(thetaA, thetaB, dAE, k, logit_pEB);

    % Use fminunc to find the optimal logit-transformed pEB that minimizes
    % KL divergence initialized at 0
    [logit_pEB, ~, ~, output] = fminunc(f, 0, options);

    % Convert logit-transformed pEB back to probability
    pEB = sigmoid(logit_pEB);

    % Warn if optimization fails and set pEB to 1
    if output.iterations == 0
        warning("Optimization failed to converge. Setting pEB = 1 as fallback.");
        pEB = 1;
    end

end

function [KL, d_logit_pEB] = pEBcost(thetaA, thetaB, dAE, k, logit_pEB)
    %pEBcost: Computes the KL divergence and its derivatives w.r.t pEB, 
    %         quantifying divergence between detection distributions with 
    %         and without Eve’s interception.
    %
    % Inputs:
    %   thetaA     - Alice’s parameters [lambdas, alpha, dAB]
    %   thetaB     - Bob’s parameters [pa0, pa1, pc0, pc1, pd0, pd1, pe]
    %   dAE        - Distance from Alice to Eve
    %   k          - Number of photons intercepted by Eve
    %   logit_pEB  - Logit-transformed transmission probability (pEB)
    %
    % Outputs:
    %   KL          - KL divergence between the two detection distributions
    %   d_logit_pEB - Gradient of KL divergence w.r.t logit-transformed pEB
    %
    % Copyright (c) 2024 Ibrahim Almosallam
    % Licensed under the MIT License (see LICENSE file for full details).

    % Convert logit-transformed pEB to probability and compute sigmoid derivative
    [pEB, ds] = sigmoid(logit_pEB);

    Delta = 1; % Assume Eve intercepts all pulses

    % Define Eve's parameters
    thetaE = {dAE, pEB, k, Delta};

    % Assume Alice and Bob use matching bases
    a = 1; % Alice's basis choice
    b = 1; % Bob's basis choice

    % Compute detection probabilities and derivatives
    [Ps1, d_thetaE1] = Pabe(thetaA, thetaB, thetaE, a, b, 1); % With interception
    Ps0 = Pabe(thetaA, thetaB, thetaE, a, b, 0);              % Without interception

    % Compute KL divergence between Ps1 and Ps0
    KL1 = Ps1 .* log(Ps1); % Expected log likelihood for Ps1
    KL0 = Ps1 .* log(Ps0); % Expected log likelihood under Ps0
    KL = sum(KL1 - KL0);   % KL divergence as sum of differences

    d_Ps1 = log(Ps1) - log(Ps0) + 1; % Derivative of KL w.r.t Ps1
    d_pEB = d_thetaE1{2};            % Derivative of Ps1 w.r.t pEB

    d_KL = sum(d_Ps1 .* d_pEB); % Derivative of KL w.r.t pEB
    d_logit_pEB = d_KL * ds;    % Derivative of KL w.r.t logit(pEB)
end