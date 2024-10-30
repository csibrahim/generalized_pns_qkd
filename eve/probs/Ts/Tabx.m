function [Ts, d_thetaE] = Tabx(thetaA, thetaB, thetaE, a, b, x)
    %Tabx: This function constructs the transition matrix for the system, 
    %      accounting for the probability of Eve intercepting a pulse, determined 
    %      by the parameter Delta. The resulting matrix reflects detection 
    %      probabilities with and without Eve's interception, and, if requested, 
    %      provides the partial derivatives with respect to the parameters in thetaE.
    %
    % Inputs:
    %   thetaA - Cell array of system parameters specific to the channel between 
    %            Alice and Bob, containing [lambdas, alpha, dAB].
    %   thetaB - Cell array of system parameters for Bob's detectors, containing 
    %            [pa0, pa1, pc0, pc1, pd0, pd1, pe].
    %   thetaE - Cell array of system parameters related to Eve’s eavesdropping, 
    %            containing [dAE, pEB, k, Delta], where Delta is the probability 
    %            that Eve intercepts the pulse.
    %   a      - Binary input representing Alice's basis choice (0 or 1).
    %   b      - Binary input representing Bob's basis choice (0 or 1).
    %   x      - Binary input representing Alice's bit choice (0 or 1).
    %
    % Outputs:
    %   Ts      - Transition matrix capturing detection probabilities for both 
    %             interception cases (Eve intercepts or does not), weighted by Delta.
    %   d_thetaE - Cell array containing the partial derivatives of the transition 
    %              matrix with respect to the parameters in thetaE.
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Retrieve the probability of Eve's interception (Delta) from thetaE
    Delta = thetaE{end};

    if nargout > 1
        % Compute transition matrices and derivatives for both interception cases
        [Ts0, d_thetaE0] = Tabxe(thetaA, thetaB, thetaE, a, b, x, 0);  % Eve does not intercept
        [Ts1, d_thetaE1] = Tabxe(thetaA, thetaB, thetaE, a, b, x, 1);  % Eve intercepts

        % Initialize cell array for derivatives
        d_thetaE = cell(1, numel(thetaE));

        % Combine the derivatives weighted by Delta for each parameter in thetaE
        for i = 1:3
            d_thetaE{i} = repmat([(1 - Delta) * d_thetaE0{i}, Delta * d_thetaE1{i}], 2, 1);
        end

        % Compute the derivative with respect to Delta
        d_thetaE{4} = repmat([-Ts0, Ts1], 2, 1);
    else
        % If no derivatives are requested, only compute the transition matrices
        Ts0 = Tabxe(thetaA, thetaB, thetaE, a, b, x, 0);
        Ts1 = Tabxe(thetaA, thetaB, thetaE, a, b, x, 1);
    end

    % Construct the overall transition matrix that accounts for both interception cases
    Ts = repmat([(1 - Delta) * Ts0, Delta * Ts1], 2, 1);
end