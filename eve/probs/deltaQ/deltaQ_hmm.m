function [deltaQs, Qs] = deltaQ_hmm(thetaA, thetaB, thetaE, both)
    %deltaQ_hmm: This function computes the probabilities of detection (Qs) and 
    %            error probabilities (deltaQs) for the hidden Markov model (HMM) case. 
    %            It considers both matching and non-matching basis configurations and 
    %            computes results across different pulse intensities.
    %
    %            Qs represents the overall detection probability. If both = true, this 
    %            includes both detectors clicking (P01 + P10 + P11). If both = false, 
    %            it represents the XOR detection probability (P01 + P10).
    %
    %            deltaQs represents the error probability, where a detection occurs at 
    %            the wrong detector.
    %
    % Inputs:
    %   thetaA - Cell array of system parameters specific to the channel between 
    %            Alice and Bob, containing [lambdas, alpha, dAB].
    %   thetaB - Cell array of system parameters for Bob's detectors, containing 
    %            [pa0, pa1, pc0, pc1, pd0, pd1, pe].
    %   thetaE - Cell array of system parameters related to Eve’s eavesdropping, 
    %            containing [dAE, pEB, k, Delta].
    %   both   - Logical flag, when true the function considers clicks from both 
    %            detectors, when false it considers XOR (default: false)
    %
    % Outputs:
    %   deltaQs - Error probabilities for each intensity and configuration
    %   Qs      - Detection probabilities for each intensity and configuration
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    if nargin < 4
        both = false;  % Default: consider XOR detection
    end

    % The emission matrix
    E = [1 0 0 0;  % p_00
         0 1 0 0;  % p_01
         0 0 1 0;  % p_10
         0 0 0 1;  % p_11
         0 1 0 0;  % p_0a
         0 0 1 0;  % p_a0
         0 0 0 1;  % p_1a
         0 0 0 1;  % p_a1
         0 0 0 1]; % p_aa

    Ns = size(E, 1);  % Number of states

    % Compute the transition matrix for the HMM
    Ts = T(thetaA, thetaB, thetaE);

    % Find the stationary distribution using the largest eigenvalue of Ts'
    [v, ~] = eigs(Ts', 1, "largestreal");
    v = abs(v);  % Ensure all values are non-negative
    sum_v = sum(v);  % Normalize the stationary distribution
    V = v / sum_v;
    
    % Reshape the stationary distribution and map it to the probabilities
    Ps = E' * reshape(V, Ns, []);

    Nl = size(Ps, 2) / (3 * 2^2);  % Number of intensities

    % Separate the probabilities for matched and non-matched bases
    not_matched = Ps(:, 1:Nl*8);

    % Merge basis probabilities for non-matching cases
    not_matched = not_matched(:, 1:4*Nl) + not_matched(:, 4*Nl+1:end);

    matched = Ps(:, (Nl*8+1):end);

    % Extract detection and error probabilities for non-matched and matched cases
    [deltaQs0, Qs0] = extract_deltaQs(not_matched, both);
    [deltaQs1, Qs1] = extract_deltaQs(matched, both);

    % Combine results for both matched and non-matched cases
    deltaQs = [deltaQs0, deltaQs1];
    Qs = [Qs0, Qs1];
end

function [deltaQs, Qs] = extract_deltaQs(M, both)
    %extract_deltaQs: This helper function computes the detection probabilities (Qs)
    %                 and error probabilities (deltaQs) for a given probability matrix M. 
    %                 It merges Alice's and Eve's bits and processes the results for the 
    %                 given 'both' flag (XOR or both detectors clicking).
    %
    % Inputs:
    %   M    - Probability matrix for the current basis configuration
    %   both - Logical flag for XOR detection or both detectors
    %
    % Outputs:
    %   deltaQs - Error probabilities for the given configuration
    %   Qs      - Detection probabilities for the given configuration

    Nl = size(M, 2) / (2^2);  % Number of intensities

    % Separate probabilities for correct and incorrect detection
    M0 = M(:, 1:2*Nl);
    M1 = M(:, 2*Nl+1:end);

    % Flip the bits for M1
    M1 = [M1(1, :); M1(3, :); M1(2, :); M1(4, :)];

    % Merge bits
    M = M0 + M1;

    % Merge Alice and Eve's bits
    M = M(:, 1:Nl) + M(:, Nl+1:end);

    % Compute the detection and error probabilities
    if both
        Qs = sum(M(2:4, :));    % Detection probabilities including both detectors
        deltaQs = sum(M(3:4, :)); % Error probabilities (wrong detector clicks)
    else
        Qs = sum(M(2:3, :));    % XOR detection probabilities
        deltaQs = M(3, :);      % Error probabilities (wrong detector clicks)
    end
end