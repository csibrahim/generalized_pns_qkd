function [Ps, d_thetaE] = P_hmm(thetaA, thetaB, thetaE)
    %P_hmm: This function computes the overall detection probabilities (Ps) 
    %       based on the stationary distribution of the hidden Markov model (HMM) 
    %       transition matrix. The stationary distribution is calculated from the 
    %       largest real eigenvalue of the transition matrix. If requested, the 
    %       function also computes the partial derivatives of the detection 
    %       probabilities with respect to the parameters in thetaE.
    %
    % Inputs:
    %   thetaA - Cell array of system parameters specific to the channel between 
    %            Alice and Bob, containing [lambdas, alpha, dAB].
    %   thetaB - Cell array of system parameters for Bob's detectors, containing 
    %            [pa0, pa1, pc0, pc1, pd0, pd1, pe].
    %   thetaE - Cell array of system parameters related to Eve’s eavesdropping, 
    %            containing [dAE, pEB, k, Delta].
    %
    % Outputs:
    %   Ps      - Array of detection probabilities [P_00, P_01, P_10, P_11, ...] 
    %             at Bob’s detectors, calculated using the stationary distribution 
    %             of the HMM.
    %   d_thetaE - Cell array of partial derivatives of the detection probabilities 
    %              with respect to the parameters in thetaE.
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    varargout = cell(1, max(nargout, 1));
    
    % E matrix maps the final state probabilities into detection probabilities
    E = [1 0 0 0;...  % p_00
         0 1 0 0;...  % p_01
         0 0 1 0;...  % p_10
         0 0 0 1;...  % p_11
         0 1 0 0;...  % p_0a
         0 0 1 0;...  % p_a0
         0 0 0 1;...  % p_1a
         0 0 0 1;...  % p_a1
         0 0 0 1];    % p_aa

    Ns = size(E, 1);  % Number of possible states

    % Retrieve the transition matrix (Ts) and optional derivatives
    [varargout{:}] = T(thetaA, thetaB, thetaE);
    Ts = varargout{1};

    % Compute the stationary distribution by finding the largest real eigenvalue
    [v, ~] = eigs(Ts', 1, "largestreal");
    v = abs(v);  % Ensure the eigenvector is non-negative
    sum_v = sum(v);
    v_norm = v / sum_v;  % Normalize the stationary distribution

    % Reshape the stationary vector to match the state space
    V = reshape(v_norm, Ns, []);
    M = E' * V;  % Map state probabilities to detection probabilities

    % Compute the detection probabilities
    Ps = foldVector(M);

    if nargout > 1
        % Identity matrix
        I = eye(size(Ts));
        
        d_thetaE = varargout{2};  % Derivatives of Ts with respect to thetaE parameters
        m = numel(d_thetaE);  % Number of parameters
        eps = 1e-15;  % Small constant for numerical stability

        % Compute the derivative of the stationary distribution w.r.t. each parameter
        for i = 1:m
            % Solve the linear system (Ts' - I) * dv/dθ = -d(Ts')/dθ * v
            dv = -(Ts' - I + eps * I) \ (d_thetaE{i}' * v);
            dv = dv / sum_v - v_norm * (sum(dv) / sum_v);  % Normalize the perturbed stationary distribution
            
            % Reshape the derivative vector and map it to detection probabilities
            dV = reshape(dv, Ns, []);
            dM = E' * dV;
            d_thetaE{i} = foldVector(dM);
        end
    end
end

% Helper function to merge and fold vectors into probabilities
function Ps = foldVector(M)
    Nl = size(M, 2) / 12;

    % Merge the unmatched cases and sum corresponding probabilities
    not_matched = M(:, 1:Nl*8);
    not_matched = not_matched(:, 1:4*Nl) + not_matched(:, 4*Nl+1:end);

    % Extract the matched case probabilities
    matched = M(:, (Nl*8+1):end);

    Ps = [extractPs(not_matched), extractPs(matched)];
end

% Helper function to extract probabilities from the matrix M
function Ps = extractPs(M)
    Nl = size(M, 2) / 4;

    % Merge the first and second halves of each row
    M = M(:, 1:2*Nl) + M(:, 2*Nl+1:end);

    % Merge Alice and Eve's probabilities
    M = M(:, 1:Nl) + M(:, Nl+1:end);

    Ps = reshape(M', 1, []);
end