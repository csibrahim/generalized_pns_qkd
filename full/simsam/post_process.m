function [thetaP, varR, thetaR, thetaR_MAP, samples, Ks, K_theory] = post_process(C, params)
    %post_process: Processes parameters, performs MAP estimation, MCMC sampling,
    %              and secure key rate calculations for QKD post-processing.
    %
    % Inputs:
    %     C         - Observed click counts
    %     params    - Cell array of parameters containing:
    %                   {true_thetas, sample_thetas, varF, varR, sigma, Nl, 
    %                    algorithm, maxIters, Ns, Nb, method, display, chunkSize}
    %
    % Outputs:
    %     thetaP       - Prior parameters for random variables
    %     varR         - Random variable names
    %     thetaR       - Random variable values (ground truth)
    %     thetaR_MAP   - MAP estimates of random variables
    %     samples      - MCMC samples from the posterior distribution of random variables
    %     Ks           - Secure key rates calculated from samples
    %     K_theory     - Theoretical secure key rates using true parameters
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % unpack parameters
    [true_thetas, sample_thetas, ...
     varF, varR, ...
     sigma, Nl, ...
     algorithm, maxIters, ...
     Ns, Nb, ...
     method, display, chunkSize] = deal(params{:});

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Process and Get Prior Parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Process fixed and random variables
    [varR, varF] = processVars(true_thetas, varF, varR, sigma);

    % Derive prior parameters and split variables into random and fixed
    [thetaP, thetaR, thetaF, varR, varF] = processParams(true_thetas, ...
                                                         sample_thetas, ...
                                                         varF, varR, sigma);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Maximum a posteriori probability (MAP) estimation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Find MAP estimation for random variables
    thetaR_MAP = MAP(C, varR, varF, thetaF, thetaP, ...
                     'algorithm', algorithm, ...
                     'maxIters', maxIters);

    % Display MAP results compared to ground truth
    printResults(thetaR_MAP, thetaR, varR);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Numerical Integration through MCMC Sampling
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Generate LaTeX-formatted labels
    param_labels = getLabels(Nl, varR);

    % Perform MCMC sampling for posterior distribution of random variables
    samples = sample(method, thetaR_MAP, Ns, Nb, ...
                     C, varR, varF, thetaF, thetaP, ...
                     'ground_truth', thetaR, ...
                     'display', display, ...
                     'labels', param_labels, ...
                     'chunkSize', chunkSize);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Estimation of Secure Key Rate and Display Results
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Calculate secure key rate, gain, and QBER using true parameters
    [EQs, Qs, Delta] = EQ(thetaR, thetaF, varR, varF);   % Error and gain probabilities
    delta = EQs ./ Qs;                                   % QBER using true parameters
    K_theory = keyRate(Delta, delta, Qs);                % Secure key rate (theoretical)

    % Calculate secure key rate, gain, and QBER using sampled parameters
    thetaRs = theta2cell(samples, varR);                 % Convert samples to cell format
    [EQs, Qs, Deltas] = EQ(thetaRs, thetaF, varR, varF); % Error and gain probabilities
    deltas = EQs ./ Qs;                                  % QBER using sampled parameters
    Ks = keyRate(Deltas, deltas, Qs);                    % Secure key rate from samples
    
    % Discard non-matching results for final key rates
    K_theory = K_theory(Nl + 1:end);
    Ks = Ks(:, Nl + 1:end);
end