function [Ps, d_thetaR] = P_iid(thetaR, thetaF, varR, varF)
    %P_iid: This function computes detection probabilities by marginalizing over 
    %       both Alice’s bit choice (x) and Eve’s interception (e) under the 
    %       assumption of independent and identically distributed (iid) pulses. 
    %       It returns detection probabilities for matching and non-matching bases 
    %       between Alice and Bob, as well as partial derivatives with respect to 
    %       the random variables specified in varR.
    %
    % Inputs:
    %   thetaR - Array of system parameters treated as random variables 
    %            (e.g., lambdas, alpha, dAB, pc0, pc1, etc.)
    %   thetaF - Array of system parameters treated as fixed variables 
    %   varR   - Cell array of variable names for the random system parameters
    %   varF   - Cell array of variable names for the fixed system parameters
    %
    % Outputs:
    %   Ps      - Detection probabilities for the cases where Alice and Bob's 
    %             bases match and where they do not match, marginalized over x and e.
    %             Ps contains two sets of probabilities: 
    %             one for the matching bases (a = 1, b = 1) and another 
    %             for non-matching bases (a = 0, b = 1).
    %   d_thetaR - Cell array of partial derivatives of the detection probabilities 
    %              with respect to the variables listed in varR.
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Check if partial derivatives are requested
    if (nargout > 1)
        % Compute probabilities and derivatives for a = 0, b = 1 (non-matching bases) 
        % and a = 1, b = 1 (matching bases)
        [Ps0, d_thetaR0] = Pab(thetaR, thetaF, varR, varF, 0, 1); % Non-matching
        [Ps1, d_thetaR1] = Pab(thetaR, thetaF, varR, varF, 1, 1); % Matching

        % Initialize d_thetaR
        d_thetaR = cell(1, numel(d_thetaR0));
    
        % Average the derivatives for a = 0, b = 1 and a = 1, b = 1
        for i = 1:numel(d_thetaR)
            d_thetaR{i} = [d_thetaR0{i}, d_thetaR1{i}] / 2;
        end

    else
        % If no derivatives are requested, just compute probabilities
        Ps0 = Pab(thetaR, thetaF, varR, varF, 0, 1); % Non-matching
        Ps1 = Pab(thetaR, thetaF, varR, varF, 1, 1); % Matching
    end

    % Combine probabilities for matching and non-matching bases
    Ps = [Ps0, Ps1] / 2;

end