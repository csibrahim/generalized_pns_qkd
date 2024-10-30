function thetaE_MAP = MAP(C, thetaA, thetaB, thetaP, varargin)
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
    %     thetaA     - Cell array of system parameters specific to the channel between 
    %                  Alice and Bob, containing [lambdas, alpha, dAB].
    %     thetaB     - Cell array of system parameters for Bob's detectors, containing 
    %                  [pa0, pa1, pc0, pc1, pd0, pd1, pe].
    %     varargin   - Optional name-value pair arguments:
    %                     'maxEpochs' - Maximum epochs for outer optimization (default: 200)
    %                     'maxIters'  - Maximum iterations for inner optimization (default: 50)
    %                     'tol'       - Tolerance threshold for early stopping (default: 1e-5)
    %
    % Outputs:
    %     thetaE_MAP - MAP estimated values for Eve's paramters
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Set up input parser with default values
    p = inputParser;
    addParameter(p, 'maxEpochs', 50);
    addParameter(p, 'maxIters', 200);
    addParameter(p, 'tol', 1e-5);

    % Parse input arguments
    parse(p, varargin{:});

    % Assign parsed values to variables
    maxEpochs = p.Results.maxEpochs;
    maxIters = p.Results.maxIters;
    tol = p.Results.tol;

    % Extract prior parameters (alphas, betas, upper bounds, lower bounds)
    [alphas, betas, ub, lb] = deal(thetaP{:});
    
    dim = 4; % number of Eve's parameters

     % Initialize Eve's paramters
    thetaE = zeros(1, dim);

    % Initialize thetaE with expected values based on the priors
    for i = 1:dim
        if any(isinf(ub(i)))
            % Unbounded variables (Gamma prior): Expected value = alpha/beta
            thetaE(i) = alphas(i) ./ betas(i) + lb(i); % Shift by lower bound
        else
            % Bounded variables (Beta prior): Expected value = alpha / (alpha + beta)
            thetaE(i) = alphas(i) ./ (alphas(i) + betas(i));
            thetaE(i) = thetaE(i) .* (ub(i) - lb(i)) + lb(i); % Scale and shift
        end
    end
    
    % Perform MAP estimation for Eve's paramters
    thetaE_MAP = MAP_E(C, thetaA, thetaB, thetaE, thetaP, maxIters, maxEpochs, tol);

end

function thetaE = MAP_E(C, thetaA, thetaB, thetaE, thetaP, maxIters, maxEpochs, tol)
    % MAP_E: Performs the MAP estimation specifically for Eve's variables by
    %        iteratively optimizing each parameter using the log-posterior PDF.
    %
    % Inputs:
    %     C          - Counts (observations) used to compute the likelihood
    %     thetaA     - Cell array of system parameters specific to the channel between 
    %                  Alice and Bob, containing [lambdas, alpha, dAB].
    %     thetaB     - Cell array of system parameters for Bob's detectors, containing 
    %                  [pa0, pa1, pc0, pc1, pd0, pd1, pe].
    %     thetaP     - Prior parameters for Eve's variables (alphas, betas, ub, lb)
    %     maxIters   - Maximum number of iterations for the inner optimization
    %     maxEpochs  - Maximum number of epochs for the outer optimization
    %     tol        - Tolerance threshold for early stopping
    %
    % Outputs:
    %     thetaE     - MAP estimated values for Eve's variables
    %
    % Copyright (c) 2024 Ibrahim Almosallam

    opt_options;  % Load default optimization options
    options.MaxIterations = maxIters;  % Set maximum iterations for fminunc

    dim = 4; % number of Eve's parameters

    for e = 1:maxEpochs
        thetaE_old = thetaE;  % Store the old value for convergence checking

        % Optimize each variable iteratively
        for i = 1:dim

            % Convert thetaE to phiE (latent variable) using the prior
            phiE = Phi(thetaE,thetaP);

            % Define the log-posterior PDF for optimization
            logpdf = @(phiE_i) logPDF_i(C, thetaA, thetaB, [phiE(1:i-1) phiE_i phiE(i+1:end)], thetaP, -1/sum(C), i); 

            % Optimize phiE_i using fminunc
            phiE_i = fminunc(logpdf, phiE(i), options);
            
            % Update the variable with the new optimized value
            phiE = [phiE(1:i-1) phiE_i phiE(i+1:end)];

            % Convert back to thetaE
            thetaE = Theta(phiE, thetaP); 

            % Display optimization progress
            progress((e-1)*dim + i, maxEpochs*dim, 'optimizing');
        end

        % Stop early if the change in thetaE is below a small tolerance
        if norm(thetaE - thetaE_old) < tol
            progress(maxEpochs*dim, maxEpochs*dim, 'optimizing');
            break;
        end
    end
end

function [log_pdf, d_pdf_i] = logPDF_i(C, thetaA, thetaB, phiE, thetaP, factor, i)

       [log_pdf, d_pdf] = logPDF(C, thetaA, thetaB, phiE, thetaP, factor);

       d_pdf_i = d_pdf(i);
end