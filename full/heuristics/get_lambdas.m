function lambdas = get_lambdas(Nl, alpha, dAB, thetaB, limit)
    % get_lambdas: Computes the optimal intensity range (lambdas) for a QKD
    %              system  to balance signal strength and maximize meaningful
    %              clicks at Bob's detectors
    %
    %              lambda_max: Maximizes XOR (single detector) clicks 
    %              lambda_min: Intensity where if Eve intercepts k = 1, she
    %                          would need pEB = 1 to remain undetected  
    %
    % Inputs:
    %   Nl     - Number of intensity levels to generate
    %   alpha  - Channel attenuation coefficient
    %   dAB    - Distance between Alice and Bob
    %   thetaB - Bob’s parameter set [pa0, pa1, pc0, pc1, pd0, pd1, pe]
    %   limit  - Upper bound for lambda intensities (default: 10)
    %
    % Outputs:
    %   lambdas - Array of linearly spaced lambda values between min and max
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Set default limit if not provided
    if nargin < 5, limit = 10; end

    % Load optimization settings
    opt_options;

    % Define the cost function for max_lambda
    f = @(logLambda) maxLambda(logLambda, alpha, dAB, thetaB);

    % Find log(lambda_max) to maximize single-click (XOR) events initialized at 0
    [log_maxLambda, ~, ~, output] = fminunc(f, 0, options);

    % Constrain to limit
    lambda_max = min(exp(log_maxLambda), limit);

    % Warn if lambda_max optimization fails and set lambda_max = limit
    if output.iterations == 0
        warning("Finding maximum lambda failed, setting lambda_max = limit = %d", limit);
        lambda_max = limit;
    end

    % Define the cost function for min_lambda
    f = @(logLambda) minLambda(logLambda, alpha, dAB, thetaB);

    % Find log(lambda_min) to ensure Eve needs perfect transmission (pEB = 1) initialized at 0
    [log_minLambda, ~, ~, output] = fminunc(f, 0, options);

    min_lambda = exp(log_minLambda);

    % Warn if min_lambda optimization fails and set lambda_min = 0
    if output.iterations == 0
        warning("Finding minimum lambda failed, setting min_lambda = %d", 0);
        min_lambda = 0;
    end

    % Ensure min_lambda <= max_lambda
    if min_lambda > lambda_max
        error("min_lambda > lambda_max ");
    end

    % Generate linearly spaced lambda values
    lambdas = linspace(min_lambda, lambda_max, Nl);
end

function [cost, d_logLambda] = maxLambda(logLambda, alpha, dAB, thetaB)
    %maxLambda: Computes the cost, no-click + double-click probabilities, 
    %           and its derivatives w.r.t. a specific intensity (lambda)
    %
    % Inputs:
    %   logLambda    - Log-transformed lambda for unconstrained optimization
    %   alpha        - Attenuation coefficient
    %   dAB          - Distance between Alice and Bob
    %   thetaB       - Bob’s parameter set [pa0, pa1, pc0, pc1, pd0, pd1, pe]
    %
    % Outputs:
    %   cost         - Cost value (no-click + double-click probabilities)
    %   d_logLambda  - Derivative of cost w.r.t. log(lambda)
    %
    % Copyright (c) 2024 Ibrahim Almosallam
    % Licensed under the MIT License (see LICENSE file for full details).


    % Define names for random and fixed variables

    % Random variables
    varR = {'dAE'};

    % Fixed variables
    varF = {'lambdas', 'alpha', 'dAB', ...
            'pa0', 'pa1', 'pc0', 'pc1', 'pd0', 'pd1', 'pe', ...
            'pEB', 'k', 'Delta'};

    lambda = exp(logLambda);
    pAB = dist2prob(dAB, alpha);

    % No Eve scenario
    dAE = 0; 
    pEB = pAB; 
    k = 0; 
    Delta = 0;

    % Define fixed parameters with Eve's settings
    thetaA = {lambda, alpha, dAB};
    thetaE = {pEB, k, Delta};
    thetaF = [thetaA thetaB thetaE];
    thetaR = {dAE};

    % Calculate probabilities and their gradients (assuming matching basis with no Eve interception)
    a = 1; % Alice'a basis
    b = 1; % Bob's basis
    e = 1; % Eve's interception flag
    [Ps, d_thetaR] = Pabe(thetaR, thetaF, varR, varF, a, b, e);

    % Cost is the sum of probabilities for no-click and double-click events
    cost = Ps(1) + Ps(end);

    dP_d_dAE = d_thetaR{1};

    % Derive dP/d_lambda from dP/d_dAE
    [~,d_pAE_d_dAE] = dist2prob(dAE,alpha);

    % dP_d_lambda = dP_d_dAE/(lambda*d_pAE_d_dAE);
    % d_logLambda = (dP_d_lambda(1) + dP_d_lambda(end)) * lambda ;

    % Similar to the above, but avoids multiplying then dividing by lambda
    d_logLambda = (dP_d_dAE(1) + dP_d_dAE(end)) / d_pAE_d_dAE ;
end

function [KL, d_logLambda] = minLambda(logLambda, alpha, dAB, thetaB)
    %minLambda: Computes KL divergence, and its derivatives w.r.t. log(lambda), 
    %           between:
    %           1) no interception.
    %           2) Eve intercepts one photon and transmits with pEB = 1.
    %
    % Inputs:
    %   logLambda    - Log-transformed lambda for optimization
    %   alpha        - Channel attenuation coefficient
    %   dAB          - Distance between Alice and Bob
    %   thetaB       - Bob’s parameter set [pa0, pa1, pc0, pc1, pd0, pd1, pe]
    %
    % Outputs:
    %   KL           - KL divergence between no-Eve and Eve interception cases
    %   d_logLambda  - Gradient of KL divergence w.r.t. log(lambda)
    %
    % Copyright (c) 2024 Ibrahim Almosallam
    % Licensed under the MIT License (see LICENSE file for full details).

    % Define names for random and fixed variables

    % Random variables
    varR = {'dAE'};

    % Fixed variables
    varF = {'lambdas', 'alpha', 'dAB', ...
            'pa0', 'pa1', 'pc0', 'pc1', 'pd0', 'pd1', 'pe', ...
            'pEB', 'k', 'Delta'};

    lambda = exp(logLambda);      % Transform to intensity
    pAB = dist2prob(dAB, alpha);  % Normal channel loss

    % Define Eve's parameters 
    dAE = 0;    % Distance between Alice and Eve
    pEB = 1;    % Perfect channel transmission
    k = 1;      % Convert log(k) back to k
    Delta = 1;  % Eve intercepts all pulses
    
    % Define parameter sets
    thetaA = {lambda, alpha, dAB};         % Alice's parameters
    thetaE = {pEB, k, Delta};              % Eve's parameters
    thetaF = [thetaA thetaB thetaE];       % Fixed parameters
    thetaR = {dAE};                        % Random parameter
    thetaS = [thetaA thetaB {pAB, 0, 0}];  % Encodes a secure channel (no Eve)

    % Assume Alice and Bob use matching bases and Eve always intercepts
    a = 1; % Alice's basis choice
    b = 1; % Bob's basis choice
    e = 1; % Eve interception flag

    [Ps1, d_thetaR1] = Pabe(thetaR, thetaF, varR, varF, a, b, e); % With Eve’s interception
    [Ps0, d_thetaR0] = Pabe(thetaR, thetaS, varR, varF, a, b, e); % Without Eve’s interception

    % Note: The effect of thetaS with e = 1 is no interception, but we do
    %       this in order to be able to extract d_lambda from d_dAE.
    %       While lambda could be defined as a random variable to retrieve 
    %       d_lambda directly, this approach aligns with the 'eve' version 
    %       for exact numerical consistency.  

    % Compute KL divergence between Ps0 and Ps1
    KL1 = Ps1 .* log(Ps1); % Expected log likelihood under Ps1
    KL0 = Ps1 .* log(Ps0); % Expected log likelihood under Ps0
    KL = sum(KL1 - KL0);   % KL divergence as the sum of differences

    % Compute the derivative of KL w.r.t dAE
    dP1_ddAE = d_thetaR1{1}; % Derivative of Ps1 w.r.t dAE
    dP0_ddAE = d_thetaR0{1}; % Derivative of Ps0 w.r.t dAE

    % Derive dP/d_lambda from dP/d_dAE
    [~,dpAE_ddAE] = dist2prob(dAE,alpha);

    % Apply chain rule to obtain the derivative of KL w.r.t lambda
    d_KL1 = (log(Ps1) - log(Ps0) + 1) .* (dP1_ddAE/dpAE_ddAE);
    d_KL0 = (Ps1 ./ Ps0) .* (dP0_ddAE/dpAE_ddAE);

    % Sum the derivative terms
    d_logLambda = sum(d_KL1 - d_KL0);
end