function [Ps, d_thetaE] = Pabx(thetaA, thetaB, thetaE, a, b, x)
    %Pabx: Computes the detection probabilities and their derivatives for all
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
    %   Ps       - Combined detection probabilities for all intensities 
    %              (lambdas) and Eve's interception flag (e)
    %   d_thetaE - Cell array of partial derivatives of the detection 
    %              probabilities with respect to the parameters in thetaE
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Extract Delta from thetaE
    Delta = thetaE{end};

    % If derivatives are requested
    if nargout > 1
        % Compute detection probabilities and derivatives for e = 0 and e = 1
        [Ps0, d_thetaE0] = Pabxe(thetaA, thetaB, thetaE, a, b, x, 0);
        [Ps1, d_thetaE1] = Pabxe(thetaA, thetaB, thetaE, a, b, x, 1);

        % Preallocate space for derivatives
        d_thetaE = cell(1, numel(thetaE));

        % Combine derivatives for each parameter weighted by Delta
        for i = 1:3
            d_thetaE{i} = Delta .* d_thetaE1{i} + (1 - Delta) .* d_thetaE0{i};
        end

        % Derivative w.r.t. Delta
        d_thetaE{end} = Ps1 - Ps0;

    else
        % Compute only the detection probabilities for e = 0 and e = 1
        Ps0 = Pabxe(thetaA, thetaB, thetaE, a, b, x, 0);
        Ps1 = Pabxe(thetaA, thetaB, thetaE, a, b, x, 1);
    end

    % Marginalize e using Delta
    Ps = Delta .* Ps1 + (1 - Delta) .* Ps0;
end
