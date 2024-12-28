function [Ps, d_thetaE] = P_hmm(thetaA, thetaB, thetaE)
    %P_hmm: Computes the detection probabilities under the HMM. assumption 
    %       and their derivatives for all intensities (lambdas), marginalized 
    %       over Eve's interception flag (e), and Alice's bit choice (x) 
    %       w.r.t. the parameters in thetaE.
    %
    % Inputs:
    %   thetaA - Alice’s parameter set [lambdas, alpha, dAB]
    %   thetaB - Bob’s parameter set [pa0, pa1, pc0, pc1, pd0, pd1, pe]
    %   thetaE - Eve’s parameter set [dAE, pEB, k, Delta]
    %
    % Outputs:
    %   Ps       - Combined detection probabilities for all intensities 
    %              (lambdas), matching and non-matching bases marginalzed 
    %              over Alice's bit choise (x) and Eve's interception flag (e)
    %   d_thetaE - Cell array of partial derivatives of detection probabilities 
    %              with respect to the parameters in thetaE
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Emission matrix
    E = sparse([1 0 0 0;...  % p_00
                0 1 0 0;...  % p_01
                0 0 1 0;...  % p_10
                0 0 0 1;...  % p_11
                0 1 0 0;...  % p_0a
                0 0 1 0;...  % p_a0
                0 0 0 1;...  % p_1a
                0 0 0 1;...  % p_a1
                0 0 0 1]);   % p_aa

    % Number of states
    Ns = size(E, 1);

    % Compute the transition matrix and its derivatives if needed
    varargout = cell(1, max(nargout, 1));
    [varargout{:}] = T(thetaA, thetaB, thetaE);

    Ts = varargout{1};  % Transition matrix
    Nt = size(Ts, 1);   % Size of the transition matrix

    % Compute the stationary distribution using the power method
    v_norm = computeStationaryDistribution(Ts, 'power');

    % Reshape the stationary vector to match the state space
    V = reshape(v_norm, Ns, []);

    % Project state probabilities onto detection probabilities
    Ps = E' * V;

    % Fold detection probabilities into matching and non-matching bases
    Ps = foldVector(Ps);

    if nargout > 1

        % Compute partial derivatives of detection probabilities w.r.t. thetaE
        d_thetaE = varargout{2};  % Transition matrix derivatives w.r.t. thetaE
        m = numel(d_thetaE);      % Number of parameters in thetaE

        % Solve the linear system (I-Ts') * dv/dθ = d(Ts')/dθ * v
        A = eye(Nt) - Ts';  
        b = cellfun(@(M) M' * v_norm, d_thetaE, 'UniformOutput', false);
        b = cell2mat(b);
        
        % Enforce the constraint sum(v) = 1
        A = [A, ones(Nt, 1); ones(1, Nt), 0];
        b = [b; zeros(1, m)];

        % Solve Ax = b using LU decomposition (faster and numerically stable)
        [L, U] = lu(A, 'vector');
        dv = U \ (L \ b);

        % Remove augmented elements introduced to enforce sum(v) = 1
        dv(end, :) = [];

        % Compute derivatives of detection probabilities w.r.t. thetaE
        for i = 1:m

            % Reshape the derivative to match the state space
            dV = reshape(dv(:, i), Ns, []);

            % Apply the same projection
            dPs = E' * dV;

            % Fold derivatives into matching and non-matching bases
            d_thetaE{i} = foldVector(dPs);

        end

    end

end

% Helper function to fold detection probabilities into matching and non-matching bases
function Ps = foldVector(Ps)

    Nl = size(Ps, 2) / 2;  % Number of intensities
    Ps0 = Ps(:, 1:Nl);     % Non-matching bases
    Ps1 = Ps(:, Nl+1:end); % Matching bases

    Ps = [reshape(Ps0', 1, []), reshape(Ps1', 1, [])];
end

% Function to compute the stationary distribution using the power or the eigenvector method
function v_norm = computeStationaryDistribution(Ts, method)
    
    if nargin < 2
        method = 'power'; % Default method if none is specified
    end

    switch method

        case 'power'  % fast for sparse matrices
            tol = 1e-10;
            max_iter = 100;
            N = size(Ts, 1);
            v_norm = ones(N, 1) / N;  % Initialize with uniform distribution

            for iter = 1:max_iter
                v_new = Ts' * v_norm;
                v_new = v_new / sum(v_new);  % Normalize
                if norm(v_new - v_norm, 1) < tol
                    break;
                end
                v_norm = v_new;
            end

        case 'eigs'  % Eigenvector method
            [v, ~] = eigs(Ts', 1, "largestreal");
            v = abs(v);           % Ensure non-negativity
            v_norm = v / sum(v);  % Normalize
    end
end
