function thetaR_MAP = MAP(C, varR, varF, thetaF, thetaP, varargin)
    %MAP: This function computes the Maximum A Posteriori (MAP) estimate for the
    %     system parameters by maximizing the posterior distribution, which is 
    %     the combination of both the likelihood (derived from observations `C`)
    %     and the prior distributions of the parameters. The function specifically
    %     optimizes Eve's parameters, while initializing the remaining variables 
    %     to the expected values of their priors and keeping them fixed during
    %     the optimization.
    %
    % Inputs:
    %     C          - Counts (observations) used to compute the likelihood
    %     varR       - Cell array of random variable names (to be optimized)
    %     varF       - Cell array of fixed variable names
    %     thetaF     - Cell array of fixed parameter values
    %     thetaP     - Cell array of prior parameters (alphas, betas, ub, lb)
    %     varargin   - Optional name-value pair arguments:
    %                     'maxEpochs' - Maximum epochs for outer optimization (default: 200)
    %                     'maxIters'  - Maximum iterations for inner optimization (default: 50)
    %                     'tol'       - Tolerance threshold for early stopping (default: 1e-5)
    %
    % Outputs:
    %     thetaR_MAP - MAP estimated values for the random variables
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Set up input parser with default values
    p = inputParser;
    addParameter(p, 'maxEpochs', 200);
    addParameter(p, 'maxIters', 50);
    addParameter(p, 'tol', 1e-5);

    % Parse input arguments
    parse(p, varargin{:});

    % Assign parsed values to variables
    maxEpochs = p.Results.maxEpochs;
    maxIters = p.Results.maxIters;
    tol = p.Results.tol;

    % Extract prior parameters (alphas, betas, upper bounds, lower bounds)
    [alphas, betas, ub, lb] = deal(thetaP{:});
    
    n = numel(alphas);  % Number of random variables
    thetaR = cell(1, n);  % Initialize cell array for random variables

    % Initialize thetaR with expected values based on the priors
    for i = 1:n
        if any(isinf(ub{i}))
            % Unbounded variables (Gamma prior): Expected value = alpha/beta
            thetaR{i} = alphas{i} ./ betas{i} + lb{i}; % Shift by lower bound
        else
            % Bounded variables (Beta prior): Expected value = alpha / (alpha + beta)
            thetaR{i} = alphas{i} ./ (alphas{i} + betas{i});
            thetaR{i} = thetaR{i} .* (ub{i} - lb{i}) + lb{i}; % Scale and shift
        end
    end
    
    % Identify and extract Eve's variables (dAE, pEB, k, Delta)
    EvesVars = {'dAE', 'pEB', 'k', 'Delta'};
    isEvesVars = ismember(varR, EvesVars);

    % Define fixed variables for the optimization (excluding Eve's variables)
    varF_E = [varF(:)', varR(~isEvesVars)];
    thetaF_E = [thetaF(:)', thetaR(~isEvesVars)];

    % Set up Eve's variables and corresponding priors
    varR_E = varR(isEvesVars);
    thetaR_E = thetaR(isEvesVars);
    alphas_E = alphas(isEvesVars);
    betas_E = betas(isEvesVars);
    ub_E = ub(isEvesVars);
    lb_E = lb(isEvesVars);
    thetaP_E = {alphas_E, betas_E, ub_E, lb_E};

    % Perform MAP estimation for Eve's variables
    thetaR_E = MAP_E(C, varR_E, varF_E, thetaR_E, thetaF_E, thetaP_E, maxIters, maxEpochs, tol);

    % Update the MAP estimate with the optimized values of Eve's variables
    thetaR(isEvesVars) = thetaR_E;
    
    % Combine all thetaR values into a single array for output
    thetaR_MAP = [thetaR{:}];
    
end

function thetaR = MAP_E(C, varR, varF, thetaR, thetaF, thetaP, maxIters, maxEpochs, tol)
    % MAP_E: Performs the MAP estimation specifically for Eve's variables by
    %        iteratively optimizing each parameter using the log-posterior PDF.
    %
    % Inputs:
    %     C          - Counts (observations) used to compute the likelihood
    %     varR       - Cell array of variable names to be optimized (Eve's variables)
    %     varF       - Cell array of fixed variable names
    %     thetaR     - Initial values for the variables to be optimized
    %     thetaF     - Values for fixed parameters during optimization
    %     thetaP     - Prior parameters for Eve's variables (alphas, betas, ub, lb)
    %     maxIters   - Maximum number of iterations for the inner optimization
    %     maxEpochs  - Maximum number of epochs for the outer optimization
    %     tol        - Tolerance threshold for early stopping
    %
    % Outputs:
    %     thetaR     - MAP estimated values for Eve's variables
    %
    % Copyright (c) 2024 Ibrahim Almosallam

    opt_options;  % Load default optimization options
    options.MaxIterations = maxIters;  % Set maximum iterations for fminunc

    % Extract prior parameters for Eve's variables
    [alphas, betas, ub, lb] = deal(thetaP{:});

    for e = 1:maxEpochs
        thetaR_old = thetaR;  % Store the old value for convergence checking

        % Optimize each variable iteratively
        for i = 1:numel(varR)
            varR_i = varR(i);  % Select current variable
            alphas_i = alphas(i);
            betas_i = betas(i);
            ub_i = ub(i);
            lb_i = lb(i);
            thetaP_i = {alphas_i, betas_i, ub_i, lb_i};

            % Define fixed variables for this iteration (exclude current variable)
            varF_i = [varF, varR(1:i-1), varR(i+1:end)];
            thetaF_i = [thetaF, thetaR(1:i-1), thetaR(i+1:end)];

            % Define the log-posterior PDF for optimization
            logpdf = @(phi) logPDF(phi, thetaF_i, varR_i, varF_i, C, thetaP_i, -1/sum(C));

            % Convert thetaR to phi (latent variable) using the prior
            phi0_i = Phi([thetaR{i}], thetaP_i);

            % Optimize phi using fminunc
            phiR_i = fminunc(logpdf, phi0_i, options);

            % Convert back to thetaR
            thetaR_i = Theta(phiR_i, thetaP_i); 

            % Update the variable with the new optimized value
            thetaR{i} = thetaR_i;

            % Display optimization progress
            progress((e-1)*numel(varR) + i, maxEpochs*numel(varR), 'optimizing');
        end

        % Stop early if the change in thetaR is below a small tolerance
        if norm([thetaR{:}] - [thetaR_old{:}]) < tol
            progress(maxEpochs*numel(varR), maxEpochs*numel(varR), 'optimizing');
            break;
        end
    end
end