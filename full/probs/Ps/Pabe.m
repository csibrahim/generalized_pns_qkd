function [Ps, d_thetaR] = Pabe(thetaR, thetaF, varR, varF, a, b, e)
    %Pabe: This function computes detection probabilities by marginalizing over x.
    %      It returns the detection probabilities and partial derivatives with 
    %      respect to the random variables specified in varR. The pulse parameters 
    %      (a, b, e) define Alice and Bob’s basis choices and whether Eve intercepted 
    %      the pulse.
    %
    % Inputs:
    %   thetaR - Array of system parameters treated as random variables 
    %            (e.g., lambdas, alpha, dAB, pc0, pc1, etc.)
    %   thetaF - Array of system parameters treated as fixed variables 
    %   varR   - Cell array of variable names for the random system parameters
    %   varF   - Cell array of variable names for the fixed system parameters
    %   a      - Binary input, Alice's basis choice for the pulse (0 or 1)
    %   b      - Binary input, Bob's basis choice for the pulse (0 or 1)
    %   e      - Binary input, Eve's interception flag 
    %            (e=1 -> she intercepted this pulse, e=0 -> she did not intercept it)
    %
    % Outputs:
    %   Ps      - Array of detection probabilities marginalized over x
    %   d_thetaR - Cell array of partial derivatives of the detection probabilities 
    %              with respect to the variables listed in varR.
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Check if partial derivatives are requested
    if nargout > 1

        % Compute probabilities and derivatives for x = 0 and x = 1
        [Ps0, d_thetaR0] = Pabxe(thetaR, thetaF, varR, varF, a, b, 0, e);
        [Ps1, d_thetaR1] = Pabxe(thetaR, thetaF, varR, varF, a, b, 1, e);
        
        % Initialize d_thetaR
        d_thetaR = cell(1, numel(d_thetaR0));
    
        % Average the derivatives for x = 0 and x = 1
        for i = 1:numel(d_thetaR)
            d_thetaR{i} = (d_thetaR0{i} + d_thetaR1{i}) / 2;
        end
        
    else

        % Compute probabilities for x = 0 and x = 1
        Ps0 = Pabxe(thetaR, thetaF, varR, varF, a, b, 0, e);
        Ps1 = Pabxe(thetaR, thetaF, varR, varF, a, b, 1, e);

    end

    % Marginalize over Alice's bit choice (x)
    Ps = (Ps0 + Ps1) / 2;

end