function [Ts, d_thetaE] = Tab(thetaA, thetaB, thetaE, a, b)
    %Tab: Computes the transition matrix and its derivatives for all
    %     intensities (lambdas), Eve's interception flag (e), and Alice's 
    %     bit choice (x) w.r.t. the parameters in thetaE. The detection 
    %     probabilities depend on Alice's and Bob's basis choices (a, b)
    %
    % Inputs:
    %   thetaA - Alice’s parameter set [lambdas, alpha, dAB]
    %   thetaB - Bob’s parameter set [pa0, pa1, pc0, pc1, pd0, pd1, pe]
    %   thetaE - Eve’s parameter set [dAE, pEB, k, Delta]
    %   a      - Binary input, Alice's basis choice (0 or 1)
    %   b      - Binary input, Bob's basis choice (0 or 1)
    %
    % Outputs:
    %   Ts       - Combined transition matrix for all intensities (lambdas),
    %              Eve's interception flag (e) and Alice's bit choice (x)
    %   d_thetaE - Cell array of partial derivatives of the transition matrix 
    %              with respect to the parameters in thetaE
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    if nargout > 1
        % Compute transition matrices and derivatives for x = 0 and x = 1
        [Ts0, d_thetaE0] = Tabx(thetaA, thetaB, thetaE, a, b, 0);
        [Ts1, d_thetaE1] = Tabx(thetaA, thetaB, thetaE, a, b, 1);

        % Preallocate for the combined derivatives
        d_thetaE = cell(1, numel(d_thetaE0));
    
        for i = 1:numel(d_thetaE)
            % Combine derivatives for each parameter (equally weighted)
            d_thetaE{i} = [d_thetaE0{i}, d_thetaE1{i}] / 2;

            % Repeat to construct the full matrix
            d_thetaE{i} = repmat(d_thetaE{i}, 2, 1);   
        end
    else
        % Compute only the transition matrices for x = 0 and x = 1
        Ts0 = Tabx(thetaA, thetaB, thetaE, a, b, 0);
        Ts1 = Tabx(thetaA, thetaB, thetaE, a, b, 1);
    end

    % Combine transition matrices (equally weighted)
    Ts = [Ts0, Ts1] / 2;

    % Repeat to construct the full matrix
    Ts = repmat(Ts, 2, 1);
end
