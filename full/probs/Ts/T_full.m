function [Ts, d_thetaR] = T_full(thetaR, thetaF, varR, varF)
    %T_full: Computes the transition matrix and its derivatives for all
    %        intensities (lambdas), Eve's interception flag (e), Alice's 
    %        bit choice (x), and Alice's and Bob's basis choices (a, b)
    %        w.r.t. the parameters in thetaR.
    %
    % Inputs:
    %   thetaR - Values of parameters listed in varR 
    %   thetaF - Values of parameters listed in varF 
    %   varR   - Cell array of random variable names
    %   varF   - Cell array of fixed variable names
    %
    % Outputs:
    %   Ts       - Combined transition matrix for all intensities (lambdas),
    %              Eve's interception flag (e), Alice's bit choice (x), and
    %              Alice's and Bob's basis choices (a, b)
    %   d_thetaR - Cell array of partial derivatives of the transition matrix 
    %              with respect to the parameters in thetaR
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).


    % Four (a, b) combinations exist: (0,0), (0,1), (1,0), (1,1).  
    % (0,0) and (1,1) share the same transition matrix, so (1,1) gets 50% weight,  
    % while (0,1) and (1,0) each get 25% weight.

    if nargout > 1
        % Compute transition matrices and derivatives for the three distinct 
        % cases: (a, b) = (0, 1), (1, 0), and (1, 1)
        [Ts01, d_thetaR01] = Tab(thetaR, thetaF, varR, varF, 0, 1);
        [Ts10, d_thetaR10] = Tab(thetaR, thetaF, varR, varF, 1, 0);
        [Ts11, d_thetaR11] = Tab(thetaR, thetaF, varR, varF, 1, 1);

        % Preallocate space for combined derivatives
        d_thetaR = cell(1, numel(d_thetaR01));
        
        for i = 1:numel(d_thetaR)
            % Combine derivatives for each parameter
            % weighted as [(0,1)/4 , (0,1)/4 , (1,1)/2] as described above
            d_thetaR{i} = [d_thetaR01{i}/4, d_thetaR10{i}/4, d_thetaR11{i}/2];

            % Repeat to construct the full matrix
            d_thetaR{i} = repmat(d_thetaR{i}, 3, 1); 
        end
    else
        % Compute transition matrices only for the three distinct 
        % cases: (a, b) = (0, 1), (1, 0), and (1, 1)
        Ts01 = Tab(thetaR, thetaF, varR, varF, 0, 1);
        Ts10 = Tab(thetaR, thetaF, varR, varF, 1, 0);
        Ts11 = Tab(thetaR, thetaF, varR, varF, 1, 1);
    end

    % Combine the transition matrices weighted as 
    % [(0,1)/4 , (0,1)/4 , (1,1)/2] as described above
    Ts = [Ts01 / 4, Ts10 / 4, Ts11 / 2];

    % Repeat to construct the full matrix
    Ts = repmat(Ts, 3, 1);

end