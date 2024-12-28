function [Ts, d_thetaR] = Tab(thetaR, thetaF, varR, varF, a, b)
    %Tab: Computes the transition matrix and its derivatives for all
    %     intensities (lambdas), Eve's interception flag (e), and Alice's 
    %     bit choice (x) w.r.t. the parameters in thetaR. The detection 
    %     probabilities depend on Alice's and Bob's basis choices (a, b)
    %
    % Inputs:
    %   thetaR - Values of parameters listed in varR 
    %   thetaF - Values of parameters listed in varF 
    %   varR   - Cell array of random variable names
    %   varF   - Cell array of fixed variable names
    %   a      - Binary input, Alice's basis choice (0 or 1)
    %   b      - Binary input, Bob's basis choice (0 or 1)
    %
    % Outputs:
    %   Ts       - Combined transition matrix for all intensities (lambdas),
    %              Eve's interception flag (e) and Alice's bit choice (x)
    %   d_thetaR - Cell array of partial derivatives of the transition matrix 
    %              with respect to the parameters in thetaR
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    if nargout > 1
        % Compute transition matrices and derivatives for x = 0 and x = 1
        [Ts0, d_thetaR0] = Tabx(thetaR, thetaF, varR, varF, a, b, 0);
        [Ts1, d_thetaR1] = Tabx(thetaR, thetaF, varR, varF, a, b, 1);

        % Preallocate for the combined derivatives
        d_thetaR = cell(1, numel(d_thetaR0));

        for i = 1:numel(d_thetaR)
            % Combine derivatives for each parameter (equally weighted)
            d_thetaR{i} = [d_thetaR0{i}, d_thetaR1{i}] / 2;

            % Repeat to construct the full matrix
            d_thetaR{i} = repmat(d_thetaR{i}, 2, 1);
        end
    else
        % Compute only the transition matrices for x = 0 and x = 1
        Ts0 = Tabx(thetaR, thetaF, varR, varF, a, b, 0);
        Ts1 = Tabx(thetaR, thetaF, varR, varF, a, b, 1);
    end

    % Combine transition matrices (equally weighted)
    Ts = [Ts0, Ts1] / 2;

    % Repeat to construct the full matrix
    Ts = repmat(Ts, 2, 1);

end