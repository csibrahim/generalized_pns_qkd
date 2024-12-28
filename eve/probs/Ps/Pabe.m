function [Ps, d_thetaE] = Pabe(thetaA, thetaB, thetaE, a, b, e)
    %Pabe: Computes detection probabilities and their derivatives for all
    %      intensities (lambdas) and marginalized over Alice's bit choice 
    %      (x) w.r.t. the parameters in thetaE. The detection probabilities 
    %      depend on Alice's and Bob's basis choices (a, b), and Eve's 
    %      interception flag (e)
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
    %   Ps       - Combined detection probabilities for all intensities 
    %              (lambdas) marginalized over Alice's bit choice (x)
    %   d_thetaE - Cell array of partial derivatives of the detection 
    %              probabilities with respect to the parameters in thetaE
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Check if partial derivatives are requested
    if nargout > 1

        % Compute detection probabilities and derivatives for x = 0 and x = 1
        [Ps0, d_thetaE0] = Pabxe(thetaA, thetaB, thetaE, a, b, 0, e);
        [Ps1, d_thetaE1] = Pabxe(thetaA, thetaB, thetaE, a, b, 1, e);

        % Preallocate space for derivatives
        d_thetaE = cell(1, numel(d_thetaE0));

        % Compute averaged partial derivatives for each parameter in thetaE
        for i = 1:numel(d_thetaE)
            d_thetaE{i} = (d_thetaE0{i} + d_thetaE1{i}) / 2;
        end

    else

        % Compute detection probabilities for x = 0 and x = 1
        Ps0 = Pabxe(thetaA, thetaB, thetaE, a, b, 0, e);
        Ps1 = Pabxe(thetaA, thetaB, thetaE, a, b, 1, e);

    end

    % Marginalize over Alice's bit choice (x)
    Ps = (Ps0 + Ps1) / 2;
end
