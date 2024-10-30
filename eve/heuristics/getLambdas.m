function lambdas = getLambdas(Nl, alpha, dAB, thetaB, limit)
    %getLambdas: Computes an optimal range of intensities (lambdas) for a QKD 
    %            system, balancing the signal strength to maximize meaningful 
    %            click events at Bob’s detectors.
    %
    %            max_lambda: The maximum intensity level, selected to maximize 
    %            informational utility by minimizing no-click and double-click 
    %            events, which contribute no useful information to Bob’s detection. 
    %            This is achieved by maximizing the XOR (single detector) click events.
    %
    %            min_lambda: The minimum intensity level, chosen such that if Eve 
    %            intercepts even a single photon, she would need perfect transmission 
    %            (pEB = 1) to remain undetectable at Bob’s detectors.
    %
    % Inputs:
    %   Nl     - Number of intensity levels to generate
    %   alpha  - Attenuation coefficient in the channel
    %   dAB    - Distance between Alice and Bob
    %   thetaB - Array of parameters associated with Bob’s detector
    %   limit  - Upper bound for lambda intensities (default: 10)
    %
    % Outputs:
    %   lambdas - Array of lambda intensities linearly spaced between min and max values
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Set default limit if not provided
    if (nargin < 5), limit = 10; end

    % Load optimization settings
    opt_options;

    % Define the cost function for max_lambda
    f = @(logLambda) maxLambda(logLambda, alpha, dAB, thetaB);

    % Use fminunc to find log(max_lambda) initialized at 0
    [log_maxLambda, ~, ~, output] = fminunc(f, 0, options);

    max_lambda = min(exp(log_maxLambda), limit); % Constrain to limit

    % Warn if optimization does not converge and set max_lambda = limit
    if output.iterations == 0
        warning("Finding maximum lambda failed, setting max_lambda = limit = %d", limit);
        max_lambda = limit;
    end

    % Define the cost function for min_lambda
    f = @(logLambda) minLambda(logLambda, alpha, dAB, thetaB);

    % Use fminunc to find log(min_lambda) initialized at 0
    [log_minLambda, ~, ~, output] = fminunc(f, 0, options);

    min_lambda = exp(log_minLambda);

    % Warn if optimization does not converge and set min_lambda = 0
    if output.iterations == 0
        warning("Finding minimum lambda failed, setting min_lambda = %d", 0);
        min_lambda = 0;
    end

    % Ensure min_lambda is not greater than max_lambda
    if min_lambda > max_lambda
        error("min_lambda > max_lambda ");
    end

    % Generate linearly spaced lambda values
    lambdas = linspace(min_lambda, max_lambda, Nl);
end

function [cost, d_logLambda] = maxLambda(logLambda, alpha, dAB, thetaB)
    %maxLambda: Computes the cost and gradient for maximizing lambda to reduce 
    %           no-click and double-click events, thereby maximizing informational 
    %           XOR clicks.
    %
    % Inputs:
    %   logLambda    - Log-transformed lambda for unconstrained optimization
    %
    % Outputs:
    %   cost         - Detection cost based on no-click and double-click events
    %   d_logLambda  - Gradient of the cost w.r.t. log(lambda)
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    lambda = exp(logLambda);
    pAB = dist2prob(dAB, alpha);

    % No Eve scenario
    dAE = 0; 
    pEB = pAB; 
    k = 0; 
    Delta = 0;

    % Define fixed parameters with Eve's settings
    thetaA = {lambda, alpha, dAB};
    thetaE = {dAE, pEB, k, Delta};
    
    % Basis and interception settings for Alice and Bob
    a = 0;
    b = 0;
    e = 1;
    [Ps, d_thetaE] = Pabe(thetaA, thetaB, thetaE, a, b, e);

    % Cost is the sum of probabilities for no-click and double-click events
    cost = Ps(1) + Ps(end);

    dP_d_dAE = d_thetaE{1};

    % Derive dP/d_lambda from dP/d_dAE
    [~,d_pAE_d_dAE] = dist2prob(dAE,alpha);

    % dP_d_lambda = dP_d_dAE/(lambda*d_pAE_d_dAE);
    % d_logLambda = (dP_d_lambda(1) + dP_d_lambda(end)) * lambda ;

    % Similar to the above, but avoids multiplying then dividing by lambda
    d_logLambda = (dP_d_dAE(1) + dP_d_dAE(end)) / d_pAE_d_dAE ;

end

function [KL, d_logLambda] = minLambda(logLambda, alpha, dAB, thetaB)
    %minLambda: Computes the KL divergence and its gradient to determine the 
    %           minimum intensity (lambda) such that if Eve intercepts at least 
    %           one photon, she would require perfect transmission (pEB = 1) to avoid 
    %           detection at Bob’s detectors.
    %
    %           The KL divergence is calculated between two detection distributions 
    %           at Bob’s detectors:
    %           1) The distribution when there is no Eve (no interception).
    %           2) The distribution when Eve intercepts one photon and transmits 
    %           with pEB = 1.
    %
    % Inputs:
    %   logLambda    - Log-transformed lambda for unconstrained optimization
    %   alpha        - Attenuation coefficient in the channel
    %   dAB          - Distance between Alice and Bob
    %   thetaB       - Array of parameters associated with Bob’s detector
    %
    % Outputs:
    %   KL           - KL divergence between the no-Eve and intercepted distributions
    %   d_logLambda  - Gradient of KL divergence w.r.t. log(lambda)
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    lambda = exp(logLambda);
    pAB = dist2prob(dAB, alpha);
    
    % Define parameters where Eve intercepts and transmits with probability pEB = 1
    dAE = 0; % Eve positioned right at Alice
    pEB = 1; % Perfect transmission
    k = 1; % Single photon interception
    Delta = 1; % Eve intercepts all pulses
    
    thetaA = {lambda, alpha, dAB};
    thetaE = {dAE, pEB, k, Delta};
    thetaS = {dAE, pAB, 0, 0}; % to encode a secure channel, no Eve
    
    % Basis choices for Alice and Bob
    a = 1;
    b = 1;

    % Compute probabilities with and without Eve's interception
    [Ps1, d_thetaE1] = Pabe(thetaA, thetaB, thetaE, a, b, 1); % With interception
    [Ps0, d_thetaE0] = Pabe(thetaA, thetaB, thetaS, a, b, 1); % Without interception

    % Calculate KL divergence between the distributions Ps1 and Ps0
    KL1 = Ps1 .* log(Ps1); % Expected log likelihood under Ps1
    KL0 = Ps1 .* log(Ps0); % Expected log likelihood under Ps0
    KL = sum(KL1 - KL0);   % KL divergence as the sum of differences

    % Compute the gradient of KL divergence w.r.t dAE
    dP1_d_dAE = d_thetaE1{1}; % Derivative of Ps1 w.r.t dAE
    dP0_d_dAE = d_thetaE0{1}; % Derivative of Ps0 w.r.t dAE

    % Derive dP/d_lambda from dP/d_dAE
    [~,d_pAE_d_dAE] = dist2prob(dAE,alpha);

    % dP1_d_lambda = dP1_d_dAE/(lambda*d_pAE_d_dAE);
    % dP0_d_lambda = dP0_d_dAE/(lambda*d_pAE_d_dAE);
    % 
    % % Apply chain rule for the gradient of KL divergence w.r.t lambda
    % d_KL1 = (log(Ps1) - log(Ps0) + 1) .* dP1_d_lambda;
    % d_KL0 = (Ps1 ./ Ps0) .* dP0_d_lambda;
    % 
    % % Sum the gradient terms and apply the lambda factor
    % d_logLambda = sum(d_KL1 - d_KL0) * lambda;

    % Apply chain rule for the gradient of KL divergence w.r.t lambda
    d_KL1 = (log(Ps1) - log(Ps0) + 1) .* (dP1_d_dAE/d_pAE_d_dAE);
    d_KL0 = (Ps1 ./ Ps0) .* (dP0_d_dAE/d_pAE_d_dAE);

    % Sum the gradient terms
    d_logLambda = sum(d_KL1 - d_KL0) ;
end