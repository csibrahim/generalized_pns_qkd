function [Ts, d_thetaE] = Tabx(thetaA, thetaB, thetaE, a, b, x)
    %Tabx: Computes the transition matrix and its derivatives for all
    %      intensities (lambdas) and Eve's interception flag (e) w.r.t. 
    %      the parameters in thetaE. The detection probabilities depend on
    %      Alice's and Bob's basis choices (a, b), Alice's bit choice (x)
    %
    % Inputs:
    %   thetaA - Alice’s parameter set [lambdas, alpha, dAB]
    %   thetaB - Bob’s parameter set [pa0, pa1, pc0, pc1, pd0, pd1, pe]
    %   thetaE - Eve’s parameter set [dAE, pEB, k, Delta]
    %   a      - Binary input, Alice's basis choice (0 or 1)
    %   b      - Binary input, Bob's basis choice (0 or 1)
    %   x      - Binary input, Alice's bit choice (0 or 1)
    %
    % Outputs:
    %   Ts       - Combined transition matrix for all intensities (lambdas)
    %              and Eve's interception flag (e)
    %   d_thetaE - Cell array of partial derivatives of the transition matrix 
    %              with respect to the parameters in thetaE
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Extract Delta (Eve's interception rate)
    Delta = thetaE{end};

    if nargout > 1
        % Compute transition matrices and derivatives for e = 0 and e = 1
        [Ts0, d_thetaE0] = Tabxe(thetaA, thetaB, thetaE, a, b, x, 0);  % No interception
        [Ts1, d_thetaE1] = Tabxe(thetaA, thetaB, thetaE, a, b, x, 1);  % Interception

        % Preallocate space for derivatives
        d_thetaE = cell(1, numel(thetaE));

        for i = 1:3
            % Combine derivatives for each parameter weighted by Delta
            d_thetaE{i} = [(1 - Delta) * d_thetaE0{i}, Delta * d_thetaE1{i}];

            % Repeat to construct the full matrix
            d_thetaE{i} = repmat(d_thetaE{i}, 2, 1);
        end

        % Derivative w.r.t. Delta
        d_thetaE{4} = [-Ts0, Ts1];

        % Repeat to construct the full matrix
        d_thetaE{4} = repmat(d_thetaE{4}, 2, 1);
    else
        % Compute only the transition matrices for e = 0 and e = 1
        Ts0 = Tabxe(thetaA, thetaB, thetaE, a, b, x, 0);
        Ts1 = Tabxe(thetaA, thetaB, thetaE, a, b, x, 1);
    end

    % Transition probabilities weighted by Delta
    Ts = [(1 - Delta) * Ts0, Delta * Ts1];

    % Repeat to construct the full matrix
    Ts = repmat(Ts, 2, 1);
end
