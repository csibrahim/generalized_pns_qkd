function thetaR_MAP = MAP(C, varR, varF, thetaF, thetaP, varargin)
    %MAP: Computes the Maximum A Posteriori (MAP) estimate for system parameters
    %     by maximizing the posterior (likelihood from `C` and priors).
    %     Optimizes random variables (`varR`), fixing others (`varF`) to `thetaF`.
    %
    % Inputs:
    %     C          - Observed click counts
    %     varR       - Cell array of random variable names
    %     varF       - Cell array of fixed variable names
    %     thetaF     - Fixed parameter values for `varF`
    %     thetaP     - Prior parameters [alphas, betas, ub, lb]
    %     varargin   - Optional name-value pairs:
    %                    'algorithm'  - trust-region or quasi-newton 
    %                                    (default: quasi-newton)
    %                    'maxIters'   - Max iterations (default: 1000)
    %
    % Outputs:
    %     thetaR_MAP - MAP estimate for the random variables
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
    algorithm = p.Results.algorithm;  % Chosen optimization algorithm
    maxIters = p.Results.maxIters;    % Max iterations for optimization

    % Extract prior parameters
    [alphas, betas, ub, lb] = deal(thetaP{:});
    
    
    dim = numel(alphas);    % Number of random variables
    thetaR = cell(1, dim);  % Initialize array for random variable estimates

    % Set initial thetaR with expected values based on the priors
    for i = 1:dim
        if any(isinf(ub{i}))
            % Semi-bounded variables (Gamma prior): E[X] = alpha / beta
            thetaR{i} = alphas{i} ./ betas{i} + lb{i}; % Shift by lower bound
        else
            % Bounded variables (Beta prior): E[X] = alpha / (alpha + beta)
            thetaR{i} = alphas{i} ./ (alphas{i} + betas{i});
            thetaR{i} = thetaR{i} .* (ub{i} - lb{i}) + lb{i}; % Scale and shift
        end
    end
    
    % Flatten thetaR
    thetaR = [thetaR{:}];

    % Define callback function for optimization progress
    callBackFunc = @(x, optimValues, state) callBack(optimValues, state, maxIters);

    % Define the objective function
    logpdf = @(phi) logPDF(phi, thetaF, varR, varF, C, thetaP, -1);
    
    % Set specified optimization options
    options.Algorithm = algorithm;                 % Optimization algorithm
    options.MaxIterations = maxIters;              % Maximum iterations
    options.MaxFunctionEvaluations = 2 * maxIters; % Max function evaluations
    options.OutputFcn = callBackFunc;              % Callback function

    % Transform initial variables from theta to phi
    phiR = Phi(thetaR, thetaP);

    % Optimize using fminunc
    phiR_MAP = fminunc(logpdf, phiR, options);

    % Transform optimal phi back to theta
    thetaR_MAP = Theta(phiR_MAP, thetaP);
end