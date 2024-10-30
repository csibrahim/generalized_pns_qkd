function max_k = getMaxK(thetaA, thetaB)
    % getMaxK: Determines the maximum number of photons, k, that Eve can intercept 
    %          from Alice's pulse such that Eve would need to use a channel 
    %          transmission probability of pEB = 1 to remain undetected.
    %
    % Inputs:
    %   thetaA - Cell array of system parameters specific to the channel between 
    %            Alice and Bob, containing [lambdas, alpha, dAB].
    %   thetaB - Cell array of system parameters for Bob's detectors, containing 
    %            [pa0, pa1, pc0, pc1, pd0, pd1, pe].
    %
    % Outputs:
    %   max_k  - Maximum value of k such that Eve requires pEB = 1 to avoid detection
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Load optimization settings
    opt_options;

    [lambdas, alpha, dAB] = thetaA{:};
    max_lambda = max(lambdas);          % Use the highest intensity in lambdas as the target intensity
    thetaA = {max_lambda, alpha, dAB};  % Redefine thetaA to only include the max intensity

    
    % Define the cost function for maximizing k
    f = @(logK) kCost(logK,thetaA, thetaB);

    % Use fminunc to find the optimal log(k) that minimizes the
    % KL divergence initialized at log(max_lambda)
    [logK,~,~,output] = fminunc(f,log(max_lambda),options);
    
    max_k = exp(logK) + 1;
    
    % Issue a warning if optimization does not converge and set max_k = max_lambda
    if output.iterations == 0
        warning("Finding maximum k failed, setting to max_k = max(max_lambda, 2)");
        max_k = max(max_lambda, 2);
    end

    % Ensure a minimum threshold for max_k
    max_k = max(max_k, 2);
    
end

function [KL, d_logK] = kCost(logK, thetaA, thetaB)
    %kCost: The KL divergence is calculated between:
    %       1) The distribution of detector clicks at Bob when Eve does not intercept.
    %       2) The distribution when Eve intercepts k photons and transmits with pEB = 1.
    %
    % Inputs:
    %   thetaA - Cell array of system parameters specific to the channel between 
    %            Alice and Bob, containing [lambdas, alpha, dAB].
    %   thetaB - Cell array of system parameters for Bob's detectors, containing 
    %            [pa0, pa1, pc0, pc1, pd0, pd1, pe].
    %
    % Outputs:
    %   KL      - KL divergence between the no-interception and interception distributions
    %   d_logK  - Gradient of the KL divergence w.r.t log(k)
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Eve's parameters
    dAE = 0;
    pEB = 1;
    k = exp(logK) + 1;
    Delta = 1;

    thetaE = {dAE, pEB, k, Delta};
    
    % Assume matching basis
    a = 1; % basis choice for Alice
    b = 1; % basis choice for Bob

    Ps0 = Pabe(thetaA, thetaB, thetaE, a, b, 0);  % No interception
    [Ps1, d_thetaE1] = Pabe(thetaA, thetaB, thetaE, a, b, 1);  % Interception with k photons

    % Compute KL divergence
    KL1 = Ps0 .* log(Ps0); % Expected log likelihood for Ps0
    KL0 = Ps0 .* log(Ps1); % Expected log likelihood under Ps1
    KL = sum(KL1 - KL0);   % Sum of differences provides KL divergence

    % Compute gradient of KL divergence w.r.t k
    d_k1 = d_thetaE1{1}; % Derivative of Ps1 w.r.t k

    % Compute gradient
    d_KL = - (Ps0 ./ Ps1) .* d_k1;
    d_logK = sum(d_KL) * (k - 1);

end