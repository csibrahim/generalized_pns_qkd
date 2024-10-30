function [Ts, d_thetaE] = T(thetaA, thetaB, thetaE)
    % T: This function computes the overall transition matrix for the system,
    %    marginalized over both Alice's and Bob's basis choices (a, b). There are
    %    four possible combinations of (a, b): 00, 01, 10, and 11. However, since 
    %    the combinations 00 and 11 produce the same transition matrix, we only 
    %    compute one matching matrix (for 11) and assign it a 50% probability. 
    %    The matrices for (a, b) = 01 and (a, b) = 10 are each assigned a 25% 
    %    probability. The function also computes the partial derivatives of the 
    %    transition matrix with respect to the parameters in thetaE, if requested.
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
    %   Ts       - Transition matrix marginalized over both Alice's and Bob's basis 
    %              choices, with appropriate probabilities for the different cases.
    %   d_thetaE - Cell array of partial derivatives of the transition matrix with 
    %              respect to the parameters in thetaE.
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    if nargout > 1
        % Compute transition matrices and derivatives for the three distinct cases:
        % (a, b) = 01, (a, b) = 10, and (a, b) = 11 (used for both 00 and 11)
        [Ts01, d_thetaE01] = Tab(thetaA, thetaB, thetaE, 0, 1);
        [Ts10, d_thetaE10] = Tab(thetaA, thetaB, thetaE, 1, 0);
        [Ts11, d_thetaE11] = Tab(thetaA, thetaB, thetaE, 1, 1);
    
        % Initialize cell array for combined derivatives
        d_thetaE = cell(1, numel(d_thetaE01));
        
        % Combine derivatives with appropriate weights for each case
        for i = 1:numel(d_thetaE)
            d_thetaE{i} = repmat([d_thetaE01{i} / 4, d_thetaE10{i} / 4, d_thetaE11{i} / 2], 3, 1); 
        end
    else
        % If no derivatives are requested, only compute the transition matrices
        Ts01 = Tab(thetaA, thetaB, thetaE, 0, 1);
        Ts10 = Tab(thetaA, thetaB, thetaE, 1, 0);
        Ts11 = Tab(thetaA, thetaB, thetaE, 1, 1);
    end

    % Combine the transition matrices, with 25% weights for Ts01 and Ts10, and 
    % 50% for Ts11 (representing both matching basis cases)
    Ts = repmat([Ts01 / 4, Ts10 / 4, Ts11 / 2], 3, 1);
end