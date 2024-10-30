function [Ts, d_thetaE] = Tab(thetaA, thetaB, thetaE, a, b)
    %Tab: This function computes the transition matrix by marginalizing over 
    %     both Alice's bit choice, x, and the eavesdropping flag, e. If requested, 
    %     the function also computes the partial derivatives of the transition
    %     matrix with respect to the parameters in thetaE.
    %
    % Inputs:
    %   thetaA - Cell array of system parameters specific to the channel between 
    %            Alice and Bob, containing [lambdas, alpha, dAB].
    %   thetaB - Cell array of system parameters for Bob's detectors, containing 
    %            [pa0, pa1, pc0, pc1, pd0, pd1, pe].
    %   thetaE - Cell array of system parameters related to Eve’s eavesdropping, 
    %            containing [dAE, pEB, k, Delta].
    %   a      - Binary input representing Alice's basis choice (0 or 1).
    %   b      - Binary input representing Bob's basis choice (0 or 1).
    %
    % Outputs:
    %   Ts      - Transition matrix marginalized over Alice's bit choice (x), 
    %             capturing detection probabilities for both bit values equally.
    %   d_thetaE - Cell array containing the partial derivatives of the transition 
    %              matrix with respect to the parameters in thetaE.
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    if nargout > 1
        % Compute transition matrices and derivatives for x = 0 and x = 1
        [Ts0, d_thetaE0] = Tabx(thetaA, thetaB, thetaE, a, b, 0);
        [Ts1, d_thetaE1] = Tabx(thetaA, thetaB, thetaE, a, b, 1);

        % Initialize cell array for combined derivatives
        d_thetaE = cell(1, numel(d_thetaE0));
    
        % Average the derivatives for x = 0 and x = 1
        for i = 1:numel(d_thetaE)
            d_thetaE{i} = repmat([d_thetaE0{i}, d_thetaE1{i}], 2, 1) / 2;   
        end
    else
        % If no derivatives are requested, only compute the transition matrices
        Ts0 = Tabx(thetaA, thetaB, thetaE, a, b, 0);
        Ts1 = Tabx(thetaA, thetaB, thetaE, a, b, 1);
    end

    % Construct the transition matrix by averaging over x = 0 and x = 1
    Ts = repmat([Ts0, Ts1], 2, 1) / 2;
end