function [Ps, d_thetaR] = Pabx(thetaR, thetaF, varR, varF, a, b, x)
    %Pabx: This function computes detection probabilities by marginalizing over e.
    %      It returns the detection probabilities and partial derivatives with 
    %      respect to the random variables specified in varR. The pulse parameters 
    %      (a, b, x) define the basis choices of Alice and Bob and Alice's bit choice.
    %
    % Inputs:
    %   thetaR - Array of system parameters treated as random variables 
    %            (e.g., lambdas, alpha, dAB, pc0, pc1, etc.)
    %   thetaF - Array of system parameters treated as fixed variables 
    %   varR   - Cell array of variable names for the random system parameters
    %   varF   - Cell array of variable names for the fixed system parameters
    %   a      - Binary input, Alice's basis choice for the pulse (0 or 1)
    %   b      - Binary input, Bob's basis choice for the pulse (0 or 1)
    %   x      - Binary input, Alice's bit choice for the pulse (0 or 1)
    %
    % Outputs:
    %   Ps      - Array of detection probabilities marginalized over e
    %   d_thetaR - Cell array of partial derivatives of the detection probabilities 
    %              with respect to the variables listed in varR.
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Combine random and fixed variables
    vars = [varR(:)' varF(:)'];
    theta = [thetaR(:)' thetaF(:)'];

    % Find the index of 'Delta' in the parameter list
    [~, ind_Delta] = ismember('Delta', vars);
    Delta = theta{ind_Delta};  % Retrieve Delta value

    % Check if partial derivatives are requested
    if nargout > 1
        % Compute probabilities and derivatives for e = 0 and e = 1
        [Ps0, d_thetaR0] = Pabxe(thetaR, thetaF, varR, varF, a, b, x, 0);
        [Ps1, d_thetaR1] = Pabxe(thetaR, thetaF, varR, varF, a, b, x, 1);

        % Initialize d_thetaR
        d_thetaR = cell(1, numel(d_thetaR0));
        
        % Combine the derivatives
        for i = 1:numel(d_thetaR)
            d_thetaR{i} = Delta .* d_thetaR1{i} + (1 - Delta) .* d_thetaR0{i};
        end

        % Handle the case where 'Delta' is in varR
        [~, ind_Delta_in_varR] = ismember('Delta', varR);
        
        if ind_Delta_in_varR
            [in_varR, ind_lambda] = ismember('lambdas', varR);

            % Adjust index if 'lambdas' is present before 'Delta'
            if in_varR && ind_lambda < ind_Delta_in_varR
                lambdas = thetaR{ind_lambda};
                Nl = size(lambdas, 2);
                ind_Delta_in_varR = ind_Delta_in_varR + Nl - 1;
            end

            % Update the derivative with respect to 'Delta'
            d_thetaR{ind_Delta_in_varR} = Ps1 - Ps0;
        end

    else
        % If no derivatives are requested, just compute probabilities
        Ps0 = Pabxe(thetaR, thetaF, varR, varF, a, b, x, 0);
        Ps1 = Pabxe(thetaR, thetaF, varR, varF, a, b, x, 1);
    end

    % Marginalize over e using Delta
    Ps = Delta .* Ps1 + (1 - Delta) .* Ps0;

end