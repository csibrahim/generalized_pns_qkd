function [Ps, d_thetaE] = P_iid(thetaA, thetaB, thetaE)
    %P_iid: Computes the detection probabilities under the i.i.d. assumption 
    %       and their derivatives for all intensities (lambdas), marginalized 
    %       over Eve's interception flag (e), and Alice's bit choice (x) 
    %       w.r.t. the parameters in thetaE.
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

    % Check if partial derivatives are requested
    if nargout > 1
        % Compute detection probabilities and derivatives for 
        % non-matching (a = 0, b = 1) and matching bases (a = 1, b = 1)
        [Ps0, d_thetaE0] = Pab(thetaA, thetaB, thetaE, 0, 1); % Non-matching basis
        [Ps1, d_thetaE1] = Pab(thetaA, thetaB, thetaE, 1, 1); % Matching basis

        % Preallocate space for combined derivatives
        d_thetaE = cell(1, numel(d_thetaE0));

        % Average the derivatives for for matching and non-matching bases
        for i = 1:numel(d_thetaE)
            d_thetaE{i} = [d_thetaE0{i}, d_thetaE1{i}] / 2;
        end
    else
        % Compute detection probabilities for matching and non-matching bases
        Ps0 = Pab(thetaA, thetaB, thetaE, 0, 1); % Non-matching basis
        Ps1 = Pab(thetaA, thetaB, thetaE, 1, 1); % Matching basis
    end

    % Combine probabilities for matching and non-matching bases and normalize
    Ps = [Ps0, Ps1] / 2;

end
