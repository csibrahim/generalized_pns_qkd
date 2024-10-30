function lambdas = getLambdas(Nl, alpha, dAB, thetaB, limit)
    % getLambdas: Computes an optimal range of intensities (lambdas) for a QKD system,
    %             balancing the signal strength to maximize meaningful click events at Bob’s
    %             detectors.
    %
    %   max_lambda: The maximum intensity level, selected to maximize informational utility
    %               by minimizing no-click and double-click events, which contribute no 
    %               useful information to Bob’s detection. This is achieved by maximizing 
    %               the XOR (single detector) click events.
    %
    %   min_lambda: The minimum intensity level, chosen such that if Eve intercepts even 
    %               a single photon, she would need perfect transmission (pEB = 1) to remain 
    %               undetectable at Bob’s detectors.
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
    if nargin < 5, limit = 10; end

    % Define variable names for random and fixed parameters
    varR = {'lambdas'}; % Random variables
    varF = {'alpha', 'dAB', ...
            'pa0', 'pa1', 'pc0', 'pc1', 'pd0', 'pd1', 'pe', ...
            'dAE', 'pEB', 'k', 'Delta'}; % Fixed variables

    % Load optimization settings
    opt_options;

    % Define the cost function for max_lambda
    f = @(logLambda) maxLambda(logLambda, varR, varF, alpha, dAB, thetaB);

    % Use fminunc to find log(max_lambda) initialized at 0
    [log_maxLambda, ~, ~, output] = fminunc(f, 0, options);

    max_lambda = min(exp(log_maxLambda), limit); % Constrain to limit

    % Warn if optimization does not converge and set max_lambda = limit
    if output.iterations == 0
        warning("Finding maximum lambda failed, setting max_lambda = limit = %d", limit);
        max_lambda = limit;
    end

    % Define the cost function for min_lambda
    f = @(logLambda) minLambda(logLambda, varR, varF, alpha, dAB, thetaB);

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

function [cost, d_logLambda] = maxLambda(logLambda, varR, varF, alpha, dAB, thetaB)
    % MAXLAMBDA: Computes the cost and gradient for maximizing lambda to reduce no-click 
    %            and double-click events, thereby maximizing informational XOR clicks.
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

    % No Eve scenario
    dAE = dAB; 
    pEB = 1; 
    k = 0; 
    Delta = 0;

    % Define fixed parameters with Eve's settings
    thetaE = {dAE, pEB, k, Delta};
    thetaF = [{alpha} {dAB} thetaB thetaE];
    thetaR = {lambda};

    % Calculate probabilities and their gradients (assuming matching basis with no Eve interception)
    a = 0; % Alice'a basis
    b = 0; % Bob's basis
    e = 1; % Eve's interception flag
    [Ps, d_thetaR] = Pabe(thetaR, thetaF, varR, varF, a, b, e);

    % Cost is the sum of probabilities for no-click and full-click events
    cost = Ps(1) + Ps(end);

    % Gradient w.r.t lambda
    d_lambda = d_thetaR{1};
    d_logLambda = (d_lambda(1) + d_lambda(end)) * lambda;
end

function [KL, d_logLambda] = minLambda(logLambda, varR, varF, alpha, dAB, thetaB)
    % MINLAMBDA: Computes the KL divergence and its gradient to determine the 
    %            minimum intensity (lambda) such that if Eve intercepts at least 
    %            one photon, she would require perfect transmission (pEB = 1) to avoid 
    %            detection at Bob’s detectors.
    %
    %            The KL divergence is calculated between two detection distributions 
    %            at Bob’s detectors:
    %            1) The distribution when there is no Eve (no interception).
    %            2) The distribution when Eve intercepts one photon and transmits with pEB = 1.
    %
    % Inputs:
    %   logLambda    - Log-transformed lambda for unconstrained optimization
    %
    % Outputs:
    %   KL           - KL divergence between the no-Eve and intercepted distributions
    %   d_logLambda  - Gradient of KL divergence w.r.t. log(lambda)
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    lambda = exp(logLambda);

    % Define parameters where Eve intercepts and transmits with probability pEB = 1
    dAE = 0;
    pEB = 1;
    k = 1;
    Delta = 1;
    
    thetaE = {dAE, pEB, k, Delta};
    thetaF = [{alpha} {dAB} thetaB thetaE];
    thetaR = {lambda};

    % Calculate probabilities and gradients for both intercepted and
    % non-intercepted cases (assuming matching basis a = 1 and b = 1)
    a = 1; % basis choice for Alice
    b = 1; % basis choice for Bob

    [Ps1, d_thetaR1] = Pabe(thetaR, thetaF, varR, varF, a, b, 1); % With Eve’s interception
    [Ps0, d_thetaR0] = Pabe(thetaR, thetaF, varR, varF, a, b, 0); % Without Eve’s interception

    % Compute KL divergence between the two probability distributions Ps1 and Ps0
    KL1 = Ps1 .* log(Ps1); % Expected log likelihood for Ps1
    KL0 = Ps1 .* log(Ps0); % Expected log likelihood under Ps0
    KL = sum(KL1 - KL0);   % Sum of differences provides KL divergence

    % Compute gradient of KL divergence w.r.t lambda
    d_lambda1 = d_thetaR1{1}; % Derivative of Ps1 w.r.t lambda
    d_lambda0 = d_thetaR0{1}; % Derivative of Ps0 w.r.t lambda

    % Apply chain rule for the gradient of KL divergence w.r.t lambda
    d_KL1 = (log(Ps1) - log(Ps0) + 1) .* d_lambda1;
    d_KL0 = (Ps1 ./ Ps0) .* d_lambda0;

    % Sum the gradient terms and apply the lambda derivative factor
    d_logLambda = sum(d_KL1 - d_KL0) * lambda;
end