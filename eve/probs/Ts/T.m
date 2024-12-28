function [Ts, d_thetaE] = T(thetaA, thetaB, thetaE)
    %T: Computes the transition matrix and its derivatives for both
    %   matching and non-matching bases (m) w.r.t. the parameters in
    %   thetaE marginalzing over Alice's bit choise (x) and Eve's
    %   interception flag (e).
    %
    % Inputs:
    %   thetaA - Alice’s parameter set [lambdas, alpha, dAB]
    %   thetaB - Bob’s parameter set [pa0, pa1, pc0, pc1, pd0, pd1, pe]
    %   thetaE - Eve’s parameter set [dAE, pEB, k, Delta]
    %
    % Outputs:
    %   Ts       - Combined transition matrix for matching and 
    %              non-matching bases (m) marginalzed over Alice's bit 
    %              choise (x) and Eve's interception flag (e)
    %   d_thetaE - Cell array of partial derivatives of the transition matrix 
    %              with respect to the parameters in thetaE
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    if nargout > 1
        % Compute transition matrices and derivatives for matching and 
        % non-matching bases
        [Ts0, d_thetaE0] = Tm(thetaA, thetaB, thetaE, 0); % non-matching
        [Ts1, d_thetaE1] = Tm(thetaA, thetaB, thetaE, 1); % matching
    
        % Preallocate space for combined derivatives
        d_thetaE = cell(1, numel(d_thetaE0));
        
        for i = 1:numel(d_thetaE)
            % Combine derivatives for each parameter (equally weighted)
            d_thetaE{i} = [d_thetaE0{i}, d_thetaE1{i}] / 2;

            % Repeat to construct the full matrix
            d_thetaE{i} = repmat(d_thetaE{i}, 2, 1); 
        end
    else
        % Compute transition matrices for matching and non-matching bases
        Ts0 = Tm(thetaA, thetaB, thetaE, 0);
        Ts1 = Tm(thetaA, thetaB, thetaE, 1);
    end

    % Combine transition matrices (equally weighted)
    Ts = [Ts0, Ts1] / 2;

    % Repeat to construct the full matrix
    Ts = repmat(Ts, 2, 1);
end
