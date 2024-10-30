function [Ps, d_thetaE] = Pabx(thetaA, thetaB, thetaE, a, b, x)
    %Pabx: This function computes the detection probabilities at Bob's detectors 
    %      by marginalizing over the eavesdropping flag, e. It returns the detection 
    %      probabilities for Bob's two detectors (00, 01, 10, 11) and, if requested, 
    %      the partial derivatives with respect to parameters in thetaE.
    %
    % Inputs:
    %   thetaA - Cell array of system parameters specific to the channel between 
    %            Alice and Bob, containing [lambdas, alpha, dAB].
    %   thetaB - Cell array of system parameters for Bob's detectors, containing 
    %            [pa0, pa1, pc0, pc1, pd0, pd1, pe].
    %   thetaE - Cell array of parameters related to Eve’s eavesdropping, 
    %            containing [dAE, pEB, k, Delta], where Delta is the weight parameter 
    %            for marginalizing over e.
    %   a      - Binary input, Alice's basis choice for the pulse (0 or 1).
    %   b      - Binary input, Bob's basis choice for the pulse (0 or 1).
    %   x      - Binary input, Alice's bit choice for the pulse (0 or 1).
    %
    % Outputs:
    %   Ps        - Array of detection probabilities [P_00, P_01, P_10, P_11] 
    %               at Bob's detectors, marginalized over e.
    %   d_thetaE  - Cell array of partial derivatives of the detection 
    %               probabilities with respect to the parameters in thetaE.
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Extract Delta from thetaE for marginalization
    Delta = thetaE{end};

    % If partial derivatives are requested
    if nargout > 1
        
        % Compute probabilities and derivatives for e = 0 and e = 1
        [Ps0, d_thetaE0] = Pabxe(thetaA, thetaB, thetaE, a, b, x, 0);
        [Ps1, d_thetaE1] = Pabxe(thetaA, thetaB, thetaE, a, b, x, 1);
        
        % Initialize cell array for derivatives with respect to thetaE
        d_thetaE = cell(1, numel(thetaE));
    
        % Compute combined partial derivatives for parameters in thetaE
        for i = 1:3
            d_thetaE{i} = Delta .* d_thetaE1{i} + (1 - Delta) .* d_thetaE0{i};
        end
        
        % Compute derivative with respect to Delta
        d_thetaE{end} = Ps1 - Ps0;

    else

        % Compute probabilities for e = 0 and e = 1
        Ps0 = Pabxe(thetaA, thetaB, thetaE, a, b, x, 0);
        Ps1 = Pabxe(thetaA, thetaB, thetaE, a, b, x, 1);

    end

    % Compute marginalized probabilities
    Ps = Delta .* Ps1 + (1 - Delta) .* Ps0;
end