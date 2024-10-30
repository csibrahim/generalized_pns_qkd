function [Ps, d_thetaE] = P_iid(thetaA, thetaB, thetaE)
    %P_iid: This function computes detection probabilities at Bob’s detectors 
    %       by marginalizing over both Alice’s bit choice (x) and Eve’s interception 
    %       (e) under the assumption of independent and identically distributed (iid) 
    %       pulses. It returns detection probabilities for matching and non-matching 
    %       bases between Alice and Bob, as well as partial derivatives with respect 
    %       to the parameters in thetaE if requested.
    %
    % Inputs:
    %   thetaA - Cell array of system parameters specific to the channel between 
    %            Alice and Bob, containing [lambdas, alpha, dAB].
    %   thetaB - Cell array of system parameters for Bob's detectors, containing 
    %            [pa0, pa1, pc0, pc1, pd0, pd1, pe].
    %   thetaE - Cell array of parameters related to Eve’s eavesdropping, 
    %            containing [dAE, pEB, k, Delta], where Delta is the weighting 
    %            parameter for marginalizing over e.
    %
    % Outputs:
    %   Ps       - Array of detection probabilities at Bob's detectors, marginalized 
    %              over x and e. It includes two sets of probabilities: one for the 
    %              matching bases (a = 1, b = 1) and another for non-matching bases (a = 0, b = 1).
    %   d_thetaE - Cell array of partial derivatives of the detection probabilities 
    %              with respect to the parameters in thetaE.
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % If partial derivatives are requested
    if nargout > 1
        % Compute probabilities and derivatives for a = 0, b = 1 (non-matching bases) 
        % and a = 1, b = 1 (matching bases)
        [Ps0, d_thetaE0] = Pab(thetaA, thetaB, thetaE, 0, 1); % Non-matching
        [Ps1, d_thetaE1] = Pab(thetaA, thetaB, thetaE, 1, 1); % Matching
        
        % Initialize cell array for derivatives with respect to thetaE
        d_thetaE = cell(1, numel(d_thetaE0));
        
        % Average the derivatives for non-matching and matching bases
        for i = 1:numel(d_thetaE)
            d_thetaE{i} = [d_thetaE0{i}, d_thetaE1{i}] / 2;
        end
    else
        % Compute probabilities for matching and non-matching bases
        Ps0 = Pab(thetaA, thetaB, thetaE, 0, 1); % Non-matching
        Ps1 = Pab(thetaA, thetaB, thetaE, 1, 1); % Matching
    end

    % Combine probabilities for matching and non-matching bases
    Ps = [Ps0, Ps1] / 2;
end