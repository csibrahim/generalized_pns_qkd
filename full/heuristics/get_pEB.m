function pEB = get_pEB(thetaA, thetaB, dAE, k)
    % get_pEB: Determines the optimal probability of transmission from Eve to Bob (pEB) 
    %          that minimizes the Kullback-Leibler (KL) divergence between the distributions
    %          of detection events at Bob's detectors, given two distinct settings:
    %          1. Eve intercepts `k` photons and transmits the remainder with a `pEB` probability.
    %          2. Eve does not intercept photons, allowing direct transmission.
    %
    % Inputs:
    %   thetaA - Parameters associated with Alice's transmission settings (lambdas, alpha, dAB)
    %   thetaB - Parameters associated with Bob's detection settings (pa, pc, pd, pe)
    %   dAE    - Distance from Alice to Eve
    %   k      - Number of photons Eve intercepts
    %
    % Output:
    %   pEB    - Optimal transmission probability from Eve to Bob (pEB) minimizing KL divergence
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Define variable sets
    varR = {'pEB'}; % Random variables 
    varF = {'lambdas', 'alpha', 'dAB', ...
            'pa0', 'pa1', 'pc0', 'pc1', 'pd0', 'pd1', 'pe',...
            'dAE', 'k', 'Delta'}; % Fixed variables
    
    % Set optimization options
    opt_options;
    
    % Define the cost function for KL divergence in terms of logit-transformed pEB
    f = @(logit_pEB) pEBcost(varR, varF, thetaA, thetaB, dAE, k, logit_pEB);

    % Use fminunc to find the optimal logit-transformed pEB that minimizes
    % KL divergence initialized at 0
    [logit_pEB, ~, ~, output] = fminunc(f, 0, options);

    % Convert the optimized logit value back to a probability
    pEB = sigmoid(logit_pEB);

    % Issue a warning if optimization does not converge and set pEB to 1
    if output.iterations == 0
        warning("Optimization failed to converge. Setting pEB = 1 as fallback.");
        pEB = 1;
    end

end

function [KL, d_logit_pEB] = pEBcost(varR, varF, thetaA, thetaB, dAE, k, logit_pEB)
    % PEBCOST: Computes the Kullback-Leibler (KL) divergence and its gradient with respect 
    %          to pEB, quantifying the divergence between two probability distributions of 
    %          detection events at Bob's detectors, one with Eve’s interception and one without.
    %
    % Inputs:
    %   varR       - Cell array of random variable names (only 'pEB' in this case)
    %   varF       - Cell array of fixed variable names
    %   thetaA     - Parameters associated with Alice
    %   thetaB     - Parameters associated with Bob
    %   dAE        - Distance from Alice to Eve
    %   k          - Number of photons intercepted by Eve
    %   logit_pEB  - Logit-transformed pEB (probability of Eve's transmission to Bob)
    %
    % Outputs:
    %   KL          - KL divergence between detection distributions for the current pEB
    %   d_logit_pEB - Gradient of the KL divergence w.r.t logit_pEB
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Convert logit-transformed pEB to probability and calculate sigmoid derivative
    [pEB, ds] = sigmoid(logit_pEB);

    % Define parameters for Eve’s settings
    Delta = 1;
    thetaE = {dAE, k, Delta};

    % Combine parameters, assigning them to fixed and random variables
    thetaR = {pEB};
    thetaF = [thetaA, thetaB, thetaE];

    % Calculate probabilities and gradients for detection events at Bob’s detectors
    % for the cases with and without interception (assuming matching basis a = 1 and b = 1)
    a = 1; % basis choice for Alice
    b = 1; % basis choice for Bob

    [Ps1, d_thetaR1] = Pabe(thetaR, thetaF, varR, varF, a, b, 1); % With Eve’s interception
    Ps0 = Pabe(thetaR, thetaF, varR, varF, a, b, 0); % Without Eve’s interception

    % Compute KL divergence between the two probability distributions Ps1 and Ps0
    KL1 = Ps1 .* log(Ps1); % Expected log likelihood for Ps1
    KL0 = Ps1 .* log(Ps0); % Expected log likelihood under Ps0
    KL = sum(KL1 - KL0);   % Sum of differences provides KL divergence

    % Compute gradient of KL divergence w.r.t pEB
    d_pEB1 = d_thetaR1{1}; % Derivative of Ps1 w.r.t pEB

    % Apply chain rule for the gradient of KL divergence w.r.t pEB
    d_KL1 = (log(Ps1) - log(Ps0) + 1) .* d_pEB1;
    
    % Sum the gradient terms and apply the sigmoid derivative factor
    d_logit_pEB = sum(d_KL1) * ds;
end