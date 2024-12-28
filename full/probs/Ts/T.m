function [Ts, d_thetaR] = T(thetaR, thetaF, varR, varF)
    %T: Computes the transition matrix and its derivatives for both
    %   matching and non-matching bases (m) w.r.t. the parameters in
    %   thetaR marginalzing over Alice's bit choise (x) and Eve's
    %   interception flag (e)
    %
    % Inputs:
    %   thetaR - Values of parameters listed in varR 
    %   thetaF - Values of parameters listed in varF 
    %   varR   - Cell array of random variable names
    %   varF   - Cell array of fixed variable names
    %
    % Outputs:
    %   Ts       - Combined transition matrix for matching and 
    %              non-matching bases (m) marginalzed over Alice's bit 
    %              choise (x) and Eve's interception flag (e)
    %   d_thetaR - Cell array of partial derivatives of the transition matrix 
    %              with respect to the parameters in thetaR
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    if nargout > 1
        % Compute transition matrices and derivatives for matching and 
        % non-matching bases
        [Ts0, d_thetaR0] = Tm(thetaR, thetaF, varR, varF, 0); % non-matching
        [Ts1, d_thetaR1] = Tm(thetaR, thetaF, varR, varF, 1); % matching

        % Preallocate space for combined derivatives
        d_thetaR = cell(1, numel(d_thetaR0));
        
        for i = 1:numel(d_thetaR)
            % Combine derivatives for each parameter (equally weighted)
            d_thetaR{i} = [d_thetaR0{i}, d_thetaR1{i}] / 2;

            % Repeat to construct the full matrix
            d_thetaR{i} = repmat(d_thetaR{i}, 2, 1);
        end
    else
        % Compute transition matrices for matching and non-matching bases
        Ts0 = Tm(thetaR, thetaF, varR, varF, 0);
        Ts1 = Tm(thetaR, thetaF, varR, varF, 1);
    end

    % Combine transition matrices (equally weighted)
    Ts = [Ts0, Ts1] / 2;

    % Repeat to construct the full matrix
    Ts = repmat(Ts, 2, 1);

end