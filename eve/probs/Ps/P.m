function varargout = P(thetaA, thetaB, thetaE)
    %P: This function computes the detection probabilities by selecting the 
    %   appropriate model based on the system parameters. If the after-pulse 
    %   probabilities (pa0, pa1) are greater than zero, it calls the hidden 
    %   Markov model (HMM) based function `P_hmm`. Otherwise, it assumes the 
    %   independent and identically distributed (iid) case and calls `P_iid`.
    %
    % Inputs:
    %   thetaA - Cell array of system parameters specific to the channel between 
    %            Alice and Bob, containing [lambdas, alpha, dAB].
    %   thetaB - Cell array of system parameters for Bob's detectors, containing 
    %            [pa0, pa1, pc0, pc1, pd0, pd1, pe].
    %   thetaE - Cell array of system parameters related to Eve’s eavesdropping, 
    %            containing [dAE, pEB, k, Delta].
    %
    % Outputs:
    %   varargout - Depending on the number of output arguments, this function
    %               returns detection probabilities or their partial derivatives 
    %               based on the selected model (iid or HMM).
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Initialize output cells
    varargout = cell(1, max(nargout, 1));

    % Retrieve after-pulse probabilities pa0 and pa1 from thetaB
    pa0 = thetaB{1};
    pa1 = thetaB{2};
    pa = [pa0, pa1];

    % Select the appropriate model based on the after-pulse probabilities
    if any(pa > 0)
        % Use the HMM model if after-pulse probabilities are greater than zero
        [varargout{:}] = P_hmm(thetaA, thetaB, thetaE);
    else
        % Use the iid model if there are no after-pulses
        [varargout{:}] = P_iid(thetaA, thetaB, thetaE);
    end
end