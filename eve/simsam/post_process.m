function [thetaE_MAP, samples, Ks, K_theory] = post_process(C, params)
    %post_process: Performs MAP estimation, MCMC sampling, and secure key 
    %              rate calculations for QKD post-processing.
    %
    % Inputs:
    %     C       - Observed click counts
    %     params  - Cell array of parameters containing:
    %                 {thetaA, thetaB, thetaE, thetaP, algorithm, maxIters, 
    %                  Ns, Nb, method, display, chunkSize}
    %
    % Outputs:
    %     thetaE_MAP - MAP estimate of Eve's parameters
    %     samples    - MCMC samples of Eve's parameter posterior distribution
    %     Ks         - Secure key rates calculated using sampled parameters
    %     K_theory   - Theoretical secure key rates using ground truth parameters
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Unpack parameters
    [thetaA, thetaB, thetaE, thetaP, ...
     algorithm, maxIters, ...
     Ns, Nb, method, display, chunkSize] = deal(params{:});

    lambdas = thetaA{1};  % Intensity levels
    Delta = thetaE{4};    % Fraction of intercepted pulses
    
    Nl = length(lambdas); % Number of intensity levels

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Maximum a posteriori probability (MAP) estimation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Find MAP estimation for Eve's parameters
    thetaE_MAP = MAP(C, thetaA, thetaB, thetaP, ...
                     'algorithm', algorithm, ...
                     'maxIters', maxIters);

    % Display MAP results compared to ground truth
    printResults(thetaE_MAP, thetaE);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Numerical Integration through MCMC sampling
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Generate LaTeX-formatted labels
    param_labels = {'$d_{\mathrm{AE}}$', '$p_{\mathrm{EB}}$', '$k$', '$\Delta$'};

    % Perform MCMC sampling for posterior distribution of Eve's parameters
    samples = sample(method, thetaE_MAP, Ns, Nb, ...
                     C, thetaA, thetaB, thetaP, ...
                     'ground_truth', thetaE, ...
                     'display', display, ...
                     'labels', param_labels, ...
                     'chunkSize', chunkSize);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Estimation of Secure Key Rate
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Calculate secure key rate, gain, and QBER using true parameters
    [EQs, Qs] = EQ(thetaA, thetaB, thetaE);  % Error and gain probabilities
    delta = EQs ./ Qs;                       % QBER
    K_theory = keyRate(Delta, delta, Qs);    % True secure-key rate

    % Calculate secure key rate, gain, and QBER using sampled parameters
    [EQs, Qs] = EQ(thetaA, thetaB, samples); % Error and gain probabilities
    Deltas = samples(:, end);                % Sampled interception fractions
    deltas = EQs ./ Qs;                      % QBER using samples
    Ks = keyRate(Deltas, deltas, Qs);        % Secure-key rate using samples
    
    % Discard non-matching results for final key rates
    K_theory = K_theory(Nl + 1:end);
    Ks = Ks(:, Nl + 1:end);

end
