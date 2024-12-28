function varargout = P(thetaA, thetaB, thetaE)
    %P: Computes detection probabilities and thier derivatives w.r.t. the 
    %   parameters in thetaE using the appropriate model based on the system
    %   parameters. If either after-pulse probabilities (pa0, pa1) are 
    %   greater than zero, the hidden Markov model (HMM) based function
    %   P_hmm is called. Otherwise, the independent and identically 
    %   distributed (i.i.d.) function P_iid is used.
    %
    % Inputs:
    %   thetaA - Alice’s parameter set [lambdas, alpha, dAB]
    %   thetaB - Bob’s parameter set [pa0, pa1, pc0, pc1, pd0, pd1, pe]
    %   thetaE - Eve’s parameter set [dAE, pEB, k, Delta]
    %
    % Outputs:
    %   Ps       - Combined detection probabilities for all intensities 
    %              (lambdas), matching and non-matching bases marginalzed 
    %              over Alice's bit choise (x) and Eve's interception flag (e)
    %   d_thetaE - Cell array of partial derivatives of detection probabilities 
    %              with respect to the parameters in thetaE
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Initialize output cells
    varargout = cell(1, max(nargout, 1));

    % Retrieve after-pulse probabilities from Bob’s parameter set
    pa0 = thetaB{1};  % After-pulse probability for detector D0
    pa1 = thetaB{2};  % After-pulse probability for detector D1
    pa = [pa0, pa1];  % Group after-pulse probabilities

    % Select the appropriate model based on after-pulse probabilities
    if any(pa > 0)
        % Call the HMM-based function if after-pulse probabilities are non-zero
        [varargout{:}] = P_hmm(thetaA, thetaB, thetaE);
    else
        % Call the iid-based function if there are no after-pulses
        [varargout{:}] = P_iid(thetaA, thetaB, thetaE);
    end
end
