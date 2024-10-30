function [Ps, d_thetaE] = Pabe(thetaA, thetaB, thetaE, a, b, e)
    %Pabe: This function computes the detection probabilities at Bob's detectors 
    %      by marginalizing over Alice's bit choice, x. It returns the detection 
    %      probabilities for Bob’s two detectors (00, 01, 10, 11) and, if requested, 
    %      the partial derivatives with respect to parameters in thetaE.
    %
    % Inputs:
    %   thetaA - Cell array of system parameters specific to the channel between 
    %            Alice and Bob, containing [lambdas, alpha, dAB].
    %   thetaB - Cell array of system parameters for Bob's detectors, containing 
    %            [pa0, pa1, pc0, pc1, pd0, pd1, pe].
    %   thetaE - Cell array of system parameters related to Eve’s eavesdropping, 
    %            containing [dAE, pEB, k].
    %   a      - Binary input, Alice's basis choice for the pulse (0 or 1).
    %   b      - Binary input, Bob's basis choice for the pulse (0 or 1).
    %   e      - Binary input, Eve's interception flag 
    %            (e=1 if she intercepted the pulse, e=0 otherwise).
    %
    % Outputs:
    %   Ps       - Array of detection probabilities [P_00, P_01, P_10, P_11] 
    %              at Bob's detectors, marginalized over Alice's bit choice, x.
    %   d_thetaE - Cell array of partial derivatives of the detection 
    %              probabilities with respect to the parameters in thetaE.
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % If partial derivatives are requested
    if (nargout > 1)
        
        % Compute probabilities and derivatives for x = 0 and x = 1
        [Ps0, d_thetaE0] = Pabxe(thetaA, thetaB, thetaE, a, b, 0, e);
        [Ps1, d_thetaE1] = Pabxe(thetaA, thetaB, thetaE, a, b, 1, e);
        
        % Initialize cell array for derivatives with respect to thetaE
        d_thetaE = cell(1, numel(d_thetaE0));
    
        % Compute averaged partial derivatives for parameters in thetaE
        for i = 1:numel(d_thetaE)
            d_thetaE{i} = (d_thetaE0{i} + d_thetaE1{i}) / 2;
        end

    else

        % Compute probabilities for x = 0 and x = 1
        Ps0 = Pabxe(thetaA, thetaB, thetaE, a, b, 0, e);
        Ps1 = Pabxe(thetaA, thetaB, thetaE, a, b, 1, e);

    end

    % Compute marginalized probabilities
    Ps = (Ps0 + Ps1) / 2;
end