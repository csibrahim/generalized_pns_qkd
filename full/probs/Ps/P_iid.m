function [Ps, d_thetaR] = P_iid(thetaR, thetaF, varR, varF)
    %P_iid: Computes the detection probabilities under the i.i.d. assumption 
    %       and their derivatives for all intensities (lambdas), marginalized 
    %       over Eve's interception flag (e), and Alice's bit choice (x) 
    %       w.r.t. the parameters in thetaR.
    %
    % Inputs:
    %   thetaR - Values of parameters listed in varR 
    %   thetaF - Values of parameters listed in varF 
    %   varR   - Cell array of random variable names
    %   varF   - Cell array of fixed variable names
    %
    % Outputs:
    %   Ps       - Combined detection probabilities for all intensities 
    %              (lambdas), matching and non-matching bases marginalzed 
    %              over Alice's bit choise (x) and Eve's interception flag (e)
    %   d_thetaR - Cell array of partial derivatives of detection probabilities 
    %              with respect to the parameters in thetaR
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Check if partial derivatives are requested
    if nargout > 1
        % Compute detection probabilities and derivatives for 
        % non-matching (a = 0, b = 1) and matching bases (a = 1, b = 1)
        [Ps0, d_thetaR0] = Pab(thetaR, thetaF, varR, varF, 0, 1); % Non-matching
        [Ps1, d_thetaR1] = Pab(thetaR, thetaF, varR, varF, 1, 1); % Matching

        % Preallocate space for combined derivatives
        d_thetaR = cell(1, numel(d_thetaR0));
    
        % Average the derivatives for for matching and non-matching bases
        for i = 1:numel(d_thetaR)
            d_thetaR{i} = [d_thetaR0{i}, d_thetaR1{i}] / 2;
        end
    else
        % Compute detection probabilities for matching and non-matching bases
        Ps0 = Pab(thetaR, thetaF, varR, varF, 0, 1); % Non-matching
        Ps1 = Pab(thetaR, thetaF, varR, varF, 1, 1); % Matching
    end

    % Combine probabilities for matching and non-matching bases and normalize
    Ps = [Ps0, Ps1] / 2;

end