function [deltaQs, Qs, Deltas] = deltaQ(thetaA, thetaB, thetaE, both)
    %deltaQ: This function computes the error probabilities (deltaQs) and 
    %        detection probabilities (Qs) based on whether the system follows an 
    %        independent and identically distributed (iid) model or a hidden Markov 
    %        model (HMM). The choice depends on the after-pulse probabilities (pa).
    %
    %        Additionally, the function extracts Deltas, which represent the 
    %        proportion of Eve's interceptions from the parameter set.
    %
    %        If the after-pulse probabilities (pa0, pa1) are greater than 0, the HMM 
    %        case is considered, and the function computes the probabilities based 
    %        on the hidden Markov model. Otherwise, the iid model is used.
    %
    %        The function computes these values for both matching and non-matching 
    %        basis configurations and across different pulse intensities.
    %
    % Inputs:
    %   thetaA - Cell array of system parameters specific to the channel between 
    %            Alice and Bob, containing [lambdas, alpha, dAB].
    %   thetaB - Cell array of system parameters for Bob's detectors, containing 
    %            [pa0, pa1, pc0, pc1, pd0, pd1, pe].
    %   thetaE - Cell array of system parameters related to Eve’s eavesdropping, 
    %            containing [dAE, pEB, k, Delta].
    %   both   - Logical flag, when true, the function considers clicks from both 
    %            detectors, when false it considers XOR (default: false)
    %
    % Outputs:
    %   deltaQs - Error probabilities for each intensity and configuration
    %   Qs      - Detection probabilities for each intensity and configuration
    %   Deltas  - Proportion of Eve's interceptions extracted from the parameters
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    if nargin < 5
        both = false;  % Default: consider XOR detection
    end

    % If thetaE is a matrix, convert to a cell array
    if(~iscell(thetaE))
        [n, m] = size(thetaE);
        thetaE = mat2cell(thetaE, n, ones(1, m));
    end

    
    lambdas = thetaA{1};
    pa0 = thetaB{1};
    pa1 = thetaB{2};
    
    [dAEs, pEBs, ks, Deltas] = deal(thetaE{:});

    Nl = length(lambdas);

    pa = [pa0 pa1];

    % Check if the system should follow the HMM model (non-zero after-pulse)
    if any(pa > 0)
        % HMM case

        n = length(Deltas);  % Number of samples

        % Initialize output arrays for the error and detection probabilities
        deltaQs = zeros(n, 2 * Nl);
        Qs = zeros(n, 2 * Nl);

        % Loop through each sample
        for i = 1:n
            % Construct a cell array for sample i
            thetaE_i = {dAEs(i), pEBs(i), ks(i), Deltas(i)};

            % Compute the HMM probabilities
            [deltaQs(i, :), Qs(i, :)] = deltaQ_hmm(thetaA, thetaB, thetaE_i, both);

            % Update the progress indicator
            if (n > 2)
                progress(i, n, 'Computing QBERs (delta) and Gains (Q)');
            end
        end
        
    else
        % iid case
        [deltaQs, Qs] = deltaQ_iid(thetaA, thetaB, thetaE, both);
    end
end