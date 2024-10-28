function [Ts, d_thetaR] = T(thetaR, thetaF, varR, varF)
    %T: This function computes the overall transition matrix for the system,
    %   marginalized over both Alice's and Bob's basis choices (a, b). There are
    %   four possible combinations of (a, b): 00, 01, 10, and 11. However, since 
    %   the combinations 00 and 11 produce the same transition matrix, we only 
    %   compute one matching matrix (for 11) and assign it a 50% probability. 
    %   The matrices for (a, b) = 01 and (a, b) = 10 are each assigned a 25% 
    %   probability. The function also computes the partial derivatives of the 
    %   transition matrix with respect to the random variables in varR, if requested.
    %
    % Inputs:
    %   thetaR - Array of system parameters treated as random variables (e.g., lambdas, Delta, etc.)
    %   thetaF - Array of system parameters treated as fixed variables
    %   varR   - Cell array of variable names corresponding to the random parameters
    %   varF   - Cell array of variable names corresponding to the fixed parameters
    %
    % Outputs:
    %   Ts      - Transition matrix marginalized over both Alice's and Bob's basis 
    %             choices, with appropriate probabilities for the different cases
    %   d_thetaR - (Optional) Cell array of partial derivatives of the transition 
    %              matrix with respect to the random parameters in varR
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    if nargout > 1
        % Compute the transition matrices and derivatives for the three distinct cases:
        % (a, b) = 01, (a, b) = 10, and (a, b) = 11 (used for both 00 and 11)
        [Ts01, d_thetaR01] = Tab(thetaR, thetaF, varR, varF, 0, 1);
        [Ts10, d_thetaR10] = Tab(thetaR, thetaF, varR, varF, 1, 0);
        [Ts11, d_thetaR11] = Tab(thetaR, thetaF, varR, varF, 1, 1);

        % Preallocate space for combined derivatives
        d_thetaR = cell(1, numel(d_thetaR01));
        
        % Combine the derivatives, assigning appropriate weights to each case
        for i = 1:numel(d_thetaR)
            d_thetaR{i} = repmat([d_thetaR01{i}/4, d_thetaR10{i}/4, d_thetaR11{i}/2], 3, 1); 
        end
    else
        % Compute the transition matrices for the three distinct cases
        Ts01 = Tab(thetaR, thetaF, varR, varF, 0, 1);
        Ts10 = Tab(thetaR, thetaF, varR, varF, 1, 0);
        Ts11 = Tab(thetaR, thetaF, varR, varF, 1, 1);
    end

    % Combine the transition matrices, assigning probabilities of 25% to Ts01 and Ts10,
    % and 50% to Ts11 (which represents both matching basis cases)
    Ts = repmat([Ts01/4, Ts10/4, Ts11/2], 3, 1);

end