function [Ps, d_thetaE] = Pab(thetaA, thetaB, thetaE, a, b)
    %Pab: Computes the detection probabilities and their derivatives for all
    %     intensities (lambdas), marginalized over Eve's interception 
    %     flag (e), and Alice's bit choice (x) w.r.t. the parameters in 
    %     thetaE. The detection probabilities depend on Alice's and Bob's 
    %     basis choices (a, b)
    %
    % Inputs:
    %   thetaA - Alice’s parameter set [lambdas, alpha, dAB]
    %   thetaB - Bob’s parameter set [pa0, pa1, pc0, pc1, pd0, pd1, pe]
    %   thetaE - Eve’s parameter set [dAE, pEB, k, Delta]
    %   a      - Binary input, Alice's basis choice (0 or 1)
    %   b      - Binary input, Bob's basis choice (0 or 1)
    %
    % Outputs:
    %   Ps       - Combined detection probabilities for all intensities 
    %              (lambdas), marginalized over Eve's interception flag (e) 
    %              and Alice's bit choice (x)
    %   d_thetaE - Cell array of partial derivatives of detection probabilities
    %              with respect to the parameters in thetaE
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Check if partial derivatives are requested
    if nargout > 1
        % Compute detection probabilities and derivatives for x = 0 and x = 1
        [Ps0, d_thetaE0] = Pabx(thetaA, thetaB, thetaE, a, b, 0);
        [Ps1, d_thetaE1] = Pabx(thetaA, thetaB, thetaE, a, b, 1);

        % Preallocate space for derivatives
        d_thetaE = cell(1, numel(d_thetaE0));

        % Average the derivatives for x = 0 and x = 1
        for i = 1:numel(d_thetaE)
            d_thetaE{i} = (d_thetaE0{i} + d_thetaE1{i}) / 2;
        end
    else
        % Compute detection probabilities for x = 0 and x = 1
        Ps0 = Pabx(thetaA, thetaB, thetaE, a, b, 0);
        Ps1 = Pabx(thetaA, thetaB, thetaE, a, b, 1);
    end

    % Marginalize over Alice's bit choice (x)
    Ps = (Ps0 + Ps1) / 2;
end
