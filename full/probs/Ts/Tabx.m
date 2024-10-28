function [Ts, d_thetaR] = Tabx(thetaR, thetaF, varR, varF, a, b, x)
    %Tabx: This function constructs the transition matrix that accounts for
    %      the possibility of Eve intercepting or not intercepting a pulse, 
    %      determined by the parameter Delta. The resulting matrix reflects the 
    %      detection probabilities for both cases (Eve intercepts or not), 
    %      incorporating Eve's transition between these two states. If requested, 
    %      the function also provides the partial derivatives of the transition 
    %      matrix with respect to the random variables in thetaR.
    %
    % Inputs:
    %   thetaR - Array of system parameters treated as random variables (e.g., lambdas, Delta)
    %   thetaF - Array of system parameters treated as fixed variables
    %   varR   - Cell array of variable names corresponding to the random parameters
    %   varF   - Cell array of variable names corresponding to the fixed parameters
    %   a      - Binary input representing Alice's basis choice (0 or 1)
    %   b      - Binary input representing Bob's basis choice (0 or 1)
    %   x      - Binary input representing Alice's bit choice (0 or 1)
    %
    % Outputs:
    %   Ts      - Transition matrix capturing the system dynamics with and without 
    %             Eve's interception, weighted by the interception probability Delta
    %   d_thetaR - (Optional) Cell array containing the partial derivatives of the 
    %              transition matrix with respect to the random parameters in varR
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Concatenate random and fixed variables for easier reference
    vars = [varR(:)' varF(:)'];
    theta = [thetaR(:)' thetaF(:)'];

    % Identify and extract the value of Delta (Eve's interception probability)
    [~, ind_Delta] = ismember('Delta', vars);
    Delta = theta{ind_Delta};

    if nargout > 1
        % Compute transition matrices and derivatives for cases where Eve does not intercept (e = 0)
        % and where she does intercept (e = 1)
        [Ts0, d_thetaR0] = Tabxe(thetaR, thetaF, varR, varF, a, b, x, 0);
        [Ts1, d_thetaR1] = Tabxe(thetaR, thetaF, varR, varF, a, b, x, 1);

        % Preallocate space for derivatives
        d_thetaR = cell(1, numel(d_thetaR1));
    
        % Combine the derivatives by weighting contributions from the cases e = 0 and e = 1
        for i = 1:numel(d_thetaR)
            d_thetaR{i} = repmat([(1 - Delta) * d_thetaR0{i}, Delta * d_thetaR1{i}], 2, 1);
        end

        % If Delta is one of the random variables, compute its derivative
        [~, ind_Delta] = ismember('Delta', varR);
        if ind_Delta
            [is_in_varR, ind_lambda] = ismember('lambdas', varR);

            % Adjust the index of Delta if lambdas precede it in the variable list
            if is_in_varR && ind_lambda < ind_Delta
                lambdas = thetaR{ind_lambda};
                Nl = length(lambdas);
                ind_Delta = ind_Delta + Nl - 1;
            end

            % Compute the derivative of the transition matrix w.r.t. Delta
            d_thetaR{ind_Delta} = repmat([-Ts0, Ts1], 2, 1);
        end
    else
        % If no derivatives are requested, only compute the transition matrices
        Ts0 = Tabxe(thetaR, thetaF, varR, varF, a, b, x, 0);
        Ts1 = Tabxe(thetaR, thetaF, varR, varF, a, b, x, 1);
    end

    % Construct the overall transition matrix that accounts for both interception cases
    Ts = repmat([(1 - Delta) * Ts0, Delta * Ts1], 2, 1);

end