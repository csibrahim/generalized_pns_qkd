function [Ts, d_thetaR] = Tab(thetaR, thetaF, varR, varF, a, b)
    %Tab: This function computes the transition matrix by marginalizing over 
    %     both Alice's bit choice, x, and the eavesdropping flag, e. If requested, 
    %     the function also computes the partial derivatives of the transition
    %     matrix with respect to the parameters in thetaR.
    % Inputs:
    %   thetaR - Array of system parameters treated as random variables (e.g., lambdas, Delta, etc.)
    %   thetaF - Array of system parameters treated as fixed variables
    %   varR   - Cell array of variable names corresponding to the random parameters
    %   varF   - Cell array of variable names corresponding to the fixed parameters
    %   a      - Binary input representing Alice's basis choice (0 or 1)
    %   b      - Binary input representing Bob's basis choice (0 or 1)
    %
    % Outputs:
    %   Ts      - Transition matrix marginalized over Alice's bit choice (x)
    %   d_thetaR - (Optional) Cell array of partial derivatives of the transition 
    %              matrix with respect to the random parameters in varR
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    if nargout > 1
        % Compute transition matrices and derivatives for x = 0 and x = 1
        [Ts0, d_thetaR0] = Tabx(thetaR, thetaF, varR, varF, a, b, 0);
        [Ts1, d_thetaR1] = Tabx(thetaR, thetaF, varR, varF, a, b, 1);

        % Preallocate for the combined derivatives
        d_thetaR = cell(1, numel(d_thetaR0));

        % Combine the derivatives, averaging over x = 0 and x = 1
        for i = 1:numel(d_thetaR)
            d_thetaR{i} = repmat([d_thetaR0{i}, d_thetaR1{i}], 2, 1) / 2;
        end
    else
        % Compute only the transition matrices for x = 0 and x = 1
        Ts0 = Tabx(thetaR, thetaF, varR, varF, a, b, 0);
        Ts1 = Tabx(thetaR, thetaF, varR, varF, a, b, 1);
    end

    % Construct the transition matrix by averaging over x = 0 and x = 1
    Ts = repmat([Ts0, Ts1], 2, 1) / 2;

end