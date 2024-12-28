function thetaE_MAP = MAP(C, thetaA, thetaB, thetaP, varargin)
    %MAP: Computes the Maximum A Posteriori (MAP) estimate for Eve's parameters
    %     by maximizing the posterior (likelihood from `C` and priors). Other 
    %     parameters are initialized to prior expectations and fixed.
    %
    % Inputs:
    %     C          - Observed click counts
    %     thetaA     - Alice’s parameter set [lambdas, alpha, dAB]
    %     thetaB     - Bob’s parameter set [pa0, pa1, pc0, pc1, pd0, pd1, pe].
    %     thetaP     - Priors for Eve’s parameters [alphas, betas, ub, lb]
    %     varargin   - Optional name-value pairs:
    %                    'algorithm'  - trust-region or quasi-newton 
    %                                    (default: quasi-newton)
    %                    'maxIters'   - Max iterations (default: 1000)
    %
    % Outputs:
    %     thetaE_MAP - MAP estimate for Eve's parameters
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Load default optimization options
    opt_options; 
    
    % Set up input parser with default values
    p = inputParser;
    addParameter(p, 'algorithm', 'quasi-newton');
    addParameter(p, 'maxIters', 1000);
   
    % Parse input arguments
    parse(p, varargin{:});

    % Assign parsed values to variables
    algorithm = p.Results.algorithm;
    maxIters = p.Results.maxIters;

    % Extract prior parameters
    [alphas, betas, ub, lb] = deal(thetaP{:});
    
    dim = length(alphas);   % number of Eve's parameters
    thetaE = zeros(1, dim); % Initialize Eve's paramters

    % Set initial thetaE with expected values based on the priors
    for i = 1:dim
        if any(isinf(ub(i)))
            % Semi-bounded variables (Gamma prior): E[X] = alpha / beta
            thetaE(i) = alphas(i) ./ betas(i) + lb(i); % Shift by lower bound
        else
            % Bounded variables (Beta prior): E[X] = alpha / (alpha + beta)
            thetaE(i) = alphas(i) ./ (alphas(i) + betas(i));
            thetaE(i) = thetaE(i) .* (ub(i) - lb(i)) + lb(i); % Scale and shift
        end
    end
    
    % Define callback function for optimization progress
    callBackFunc = @(x, optimValues, state) callBack(optimValues, state, maxIters);
    
    % Define the objective function
    logpdf = @(phi) logPDF(C, thetaA, thetaB, phi, thetaP, -1);
    
    % Set specified optimization options
    options.Algorithm = algorithm;                 % Optimization algorithm
    options.MaxIterations = maxIters;              % Maximum iterations
    options.MaxFunctionEvaluations = 2 * maxIters; % Max function evaluations
    options.OutputFcn = callBackFunc;              % Callback function

    % Transform initial variables from theta to phi
    phiE = Phi(thetaE, thetaP);
    
    % Optimize using fminunc
    phiE_MAP = fminunc(logpdf, phiE, options);

    % Transform optimal phi back to theta
    thetaE_MAP = Theta(phiE_MAP, thetaP);
end