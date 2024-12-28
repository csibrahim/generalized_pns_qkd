function [Ts, d_thetaR] = Tabx(thetaR, thetaF, varR, varF, a, b, x)
    %Tabx: Computes the transition matrix and its derivatives for all
    %      intensities (lambdas) and Eve's interception flag (e) w.r.t. 
    %      the parameters in thetaR. The detection probabilities depend on
    %      Alice's and Bob's basis choices (a, b), Alice's bit choice (x)
    %
    % Inputs:
    %   thetaR - Values of parameters listed in varR 
    %   thetaF - Values of parameters listed in varF 
    %   varR   - Cell array of random variable names
    %   varF   - Cell array of fixed variable names
    %   a      - Binary input, Alice's basis choice (0 or 1)
    %   b      - Binary input, Bob's basis choice (0 or 1)
    %   x      - Binary input, Alice's bit choice (0 or 1)
    %
    % Outputs:
    %   Ts       - Combined transition matrix for all intensities (lambdas)
    %              and Eve's interception flag (e)
    %   d_thetaR - Cell array of partial derivatives of the transition matrix 
    %              with respect to the parameters in thetaR
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Concatenate random and fixed variables for easier reference
    vars = [varR(:)' varF(:)'];
    theta = [thetaR(:)' thetaF(:)'];

    % Identify and extract Delta (Eve's interception rate)
    [~, ind_Delta] = ismember('Delta', vars);
    Delta = theta{ind_Delta};

    if nargout > 1
        % Compute transition matrices and derivatives for e = 0 and e = 1
        [Ts0, d_thetaR0] = Tabxe(thetaR, thetaF, varR, varF, a, b, x, 0);  % No interception
        [Ts1, d_thetaR1] = Tabxe(thetaR, thetaF, varR, varF, a, b, x, 1);  % Interception

        % Preallocate space for derivatives
        d_thetaR = cell(1, numel(d_thetaR1));
    
        for i = 1:numel(d_thetaR)
            % Combine derivatives for each parameter weighted by Delta
            d_thetaR{i} = [(1 - Delta) * d_thetaR0{i}, Delta * d_thetaR1{i}];

            % Repeat to construct the full matrix
            d_thetaR{i} = repmat(d_thetaR{i}, 2, 1);
        end

        % If Delta is one of the random variables, compute its derivative
        [~, ind_Delta] = ismember('Delta', varR);
        if ind_Delta
            [is_in_varR, ind_lambda] = ismember('lambdas', varR);

            % Adjust Delta's index if lambdas come earlier in the variable list
            if is_in_varR && ind_lambda < ind_Delta
                lambdas = thetaR{ind_lambda};
                Nl = length(lambdas);
                ind_Delta = ind_Delta + Nl - 1;
            end

            % Derivative w.r.t. Delta
            d_thetaR{ind_Delta} = [-Ts0, Ts1];

            % Repeat to construct the full matrix
            d_thetaR{ind_Delta} = repmat(d_thetaR{ind_Delta}, 2, 1);
        end
    else
        % Compute only the transition matrices for e = 0 and e = 1
        Ts0 = Tabxe(thetaR, thetaF, varR, varF, a, b, x, 0);
        Ts1 = Tabxe(thetaR, thetaF, varR, varF, a, b, x, 1);
    end

    % Transition probabilities weighted by Delta
    Ts = [(1 - Delta) * Ts0, Delta * Ts1];

    % Repeat to construct the full matrix
    Ts = repmat(Ts, 2, 1);

end