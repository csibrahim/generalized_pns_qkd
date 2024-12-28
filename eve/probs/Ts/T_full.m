function [Ts, d_thetaE] = T_full(thetaA, thetaB, thetaE)
    %T_full: Computes the transition matrix and its derivatives for all
    %        intensities (lambdas), Eve's interception flag (e), Alice's 
    %        bit choice (x), and Alice's and Bob's basis choices (a, b)
    %        w.r.t. the parameters in thetaE.
    %
    % Inputs:
    %   thetaA - Alice’s parameter set [lambdas, alpha, dAB]
    %   thetaB - Bob’s parameter set [pa0, pa1, pc0, pc1, pd0, pd1, pe]
    %   thetaE - Eve’s parameter set [dAE, pEB, k, Delta]
    %
    % Outputs:
    %   Ts       - Combined transition matrix for all intensities (lambdas),
    %              Eve's interception flag (e), Alice's bit choice (x), and
    %              Alice's and Bob's basis choices (a, b)
    %   d_thetaE - Cell array of partial derivatives of the transition matrix 
    %              with respect to the parameters in thetaE
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Four (a, b) combinations exist: (0,0), (0,1), (1,0), (1,1).  
    % (0,0) and (1,1) share the same transition matrix, so (1,1) gets 50% weight,  
    % while (0,1) and (1,0) each get 25% weight.  

    if nargout > 1
        % Compute transition matrices and derivatives for the three distinct 
        % cases: (a, b) = (0, 1), (1, 0), and (1, 1)
        [Ts01, d_thetaE01] = Tab(thetaA, thetaB, thetaE, 0, 1);
        [Ts10, d_thetaE10] = Tab(thetaA, thetaB, thetaE, 1, 0);
        [Ts11, d_thetaE11] = Tab(thetaA, thetaB, thetaE, 1, 1);
    
        % Preallocate space for combined derivatives
        d_thetaE = cell(1, numel(d_thetaE01));
        
        for i = 1:numel(d_thetaE)
            % Combine derivatives for each parameter
            % weighted as [(0,1)/4 , (0,1)/4 , (1,1)/2] as described above
            d_thetaE{i} = [d_thetaE01{i} / 4, d_thetaE10{i} / 4, d_thetaE11{i} / 2];

            % Repeat to construct the full matrix
            d_thetaE{i} = repmat(d_thetaE{i}, 3, 1); 
        end
    else
        % Compute transition matrices only for the three distinct 
        % cases: (a, b) = (0, 1), (1, 0), and (1, 1)
        Ts01 = Tab(thetaA, thetaB, thetaE, 0, 1);
        Ts10 = Tab(thetaA, thetaB, thetaE, 1, 0);
        Ts11 = Tab(thetaA, thetaB, thetaE, 1, 1);
    end

    % Combine the transition matrices weighted as 
    % [(0,1)/4 , (0,1)/4 , (1,1)/2] as described above
    Ts = [Ts01 / 4, Ts10 / 4, Ts11 / 2];

    % Repeat to construct the full matrix
    Ts = repmat(Ts, 3, 1);
end
