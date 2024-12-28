function [EQs, Qs, Deltas] = EQ(thetaA, thetaB, thetaE, both)
    %EQ: Computes signal (Q) and error (EQ) probabilities probabilities
    %    based on whether the system follows an independent and identically 
    %    distributed (i.i.d.) model or a hidden Markov model (HMM). If 
    %    either after-pulse probabilities (pa0, pa1) are greater than zero,
    %    the hidden Markov model (HMM) based function EQ_hmm is called. 
    %    Otherwise, the independent and identically distributed (i.i.d.) 
    %    function EQ_iid is used.
    %
    % Inputs:
    %   thetaA - Alice’s parameter set [lambdas, alpha, dAB]
    %   thetaB - Bob’s parameter set [pa0, pa1, pc0, pc1, pd0, pd1, pe]
    %   thetaE - Eve’s parameter set [dAE, pEB, k, Delta]
    %   both   - Logical flag (default: false):
    %            - false: Only XORs are considered signal events
    %            - true:  Includes ANDs as signal and error events
    %
    % Outputs:
    %   EQs    - Error probabilities for each intensity and configuration
    %    Qs    - Signal probabilities for each intensity and configuration
    %   Deltas - Proportion of Eve's interceptions extracted from the parameters
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    if nargin < 4
        both = false;  % Default to XOR detection
    end

    % Convert thetaE to cell array if provided as a matrix
    if ~iscell(thetaE)
        [n, m] = size(thetaE);
        thetaE = mat2cell(thetaE, n, ones(1, m));
    end

    % Retrieve after-pulse probabilities from Bob’s parameter set
    pa0 = thetaB{1};  % After-pulse probability for detector D0
    pa1 = thetaB{2};  % After-pulse probability for detector D1
    pa = [pa0, pa1];  % Group after-pulse probabilities

    % Extract Eve's parameters
    [dAEs, pEBs, ks, Deltas] = deal(thetaE{:}); 

    Nl = length(thetaA{1});  % Number of intensities

    % Determine if the HMM model should be used (non-zero after-pulse)
    if any(pa > 0) 

        % HMM case

        % Number of samples
        n = size(dAEs, 1); 

        % Initialize outputs for error and signal probabilities
        EQs = zeros(n, 2 * Nl);
         Qs = zeros(n, 2 * Nl);

        % Loop over each sample
        for i = 1:n
            % Construct a cell array for sample i
            thetaE_i = {dAEs(i), pEBs(i), ks(i), Deltas(i)};

            % Compute the HMM probabilities
            [EQs(i, :), Qs(i, :)] = EQ_hmm(thetaA, thetaB, thetaE_i, both);

            % Display progress if there are multiple samples
            if n > 2
                progress(i, n, 'Computing error and detection probabilities (Q and EQ)');
            end
        end

    else
        % iid case
        [EQs, Qs] = EQ_iid(thetaA, thetaB, thetaE, both);
    end
end
