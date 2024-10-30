function [Ps, d_thetaE] = Pab(thetaA, thetaB, thetaE, a, b)
    %Pab: This function computes the detection probabilities at Bob's detectors 
    %     by marginalizing over both Alice's bit choice, x, and the eavesdropping 
    %     flag, e. It returns the detection probabilities for Bob’s two detectors 
    %     (00, 01, 10, 11) and, if requested, the partial derivatives with respect 
    %     to parameters in thetaE.
    %
    % Inputs:
    %   thetaA - Cell array of system parameters specific to the channel between 
    %            Alice and Bob, containing [lambdas, alpha, dAB].
    %   thetaB - Cell array of system parameters for Bob's detectors, containing 
    %            [pa0, pa1, pc0, pc1, pd0, pd1, pe].
    %   thetaE - Cell array of system parameters related to Eve’s eavesdropping, 
    %            containing [dAE, pEB, k, Delta], where Delta is the weighting 
    %            parameter for marginalizing over e.
    %   a      - Binary input, Alice's basis choice for the pulse (0 or 1).
    %   b      - Binary input, Bob's basis choice for the pulse (0 or 1).
    %
    % Outputs:
    %   Ps       - Array of detection probabilities [P_00, P_01, P_10, P_11], for
    %              each intensity, at Bob's detectors, marginalized over both x and e.
    %   d_thetaE - Cell array of partial derivatives of the detection probabilities 
    %              with respect to the parameters in thetaE.
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % If partial derivatives are requested
    if nargout > 1
        
        % Compute probabilities and derivatives for x = 0 and x = 1
        [Ps0, d_thetaE0] = Pabx(thetaA, thetaB, thetaE, a, b, 0);
        [Ps1, d_thetaE1] = Pabx(thetaA, thetaB, thetaE, a, b, 1);
        
        % Initialize cell array for derivatives with respect to thetaE
        d_thetaE = cell(1, numel(d_thetaE0));
    
        % Average the derivatives for x = 0 and x = 1
        for i = 1:numel(d_thetaE)
            d_thetaE{i} = (d_thetaE0{i} + d_thetaE1{i}) / 2;
        end

    else

        % Compute probabilities for x = 0 and x = 1
        Ps0 = Pabx(thetaA, thetaB, thetaE, a, b, 0);
        Ps1 = Pabx(thetaA, thetaB, thetaE, a, b, 1);

    end

    % Compute marginalized probabilities over x
    Ps = (Ps0 + Ps1) / 2;
end