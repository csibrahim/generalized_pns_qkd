function [Ps, d_thetaR] = P_hmm(thetaR, thetaF, varR, varF)
    %P_hmm: This function computes the overall detection probabilities (Ps) 
    %       based on the stationary distribution of the hidden Markov model (HMM) 
    %       transition matrix. The stationary distribution is calculated from the 
    %       largest real eigenvalue of the transition matrix. If requested, the 
    %       function also computes the partial derivatives of the detection 
    %       probabilities with respect to the random variables in thetaR.
    %
    % Inputs:
    %   thetaR - Array of system parameters treated as random variables
    %   thetaF - Array of system parameters treated as fixed variables
    %   varR   - Cell array of variable names corresponding to the random parameters
    %   varF   - Cell array of variable names corresponding to the fixed parameters
    %
    % Outputs:
    %   Ps      - Detection probabilities [P_00, P_01, P_10, P_11, ...]
    %   d_thetaR - (Optional) Cell array of partial derivatives of the detection 
    %              probabilities with respect to the random parameters in varR
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    varargout = cell(1, nargout);

    % E matrix maps the final state probabilities into detection probabilities
    E = [1 0 0 0;... % p_00
         0 1 0 0;... % p_01
         0 0 1 0;... % p_10
         0 0 0 1;... % p_11
         0 1 0 0;... % p_0a
         0 0 1 0;... % p_a0
         0 0 0 1;... % p_1a
         0 0 0 1;... % p_a1
         0 0 0 1];   % p_aa

    Ns = size(E, 1);  % Number of possible states

    % Retrieve the transition matrix (Ts) and optional derivatives
    [varargout{:}] = T(thetaR, thetaF, varR, varF);
    Ts = varargout{1};

    % Compute the stationary distribution by finding the largest real eigenvalue
    [v, ~] = eigs(Ts', 1, "largestreal");
    v = abs(v);  % Ensure the eigenvector is non-negative
    v_norm = v / sum(v);  % Normalize the stationary distribution

    % Reshape the stationary vector to match the state space
    V = reshape(v_norm, Ns, []);
    M = E' * V;  % Map state probabilities to detection probabilities

    % Compute the detection probabilities
    Ps = foldVector(M);

    if nargout > 1
        
        % Compute derivatives with respect to the random variables in thetaR
        I = eye(size(Ts));  % Identity matrix
        d_thetaR = varargout{2};
        m = numel(d_thetaR);
        eps = 1e-15;  % Small constant for numerical stability

        % Compute the derivative of the stationary distribution w.r.t. each parameter
        for i = 1:m
            % Solve the linear system (Ts' - I) * dv/dθ = -d(Ts')/dθ * v
            dv = -(Ts' - I + eps * I) \ (d_thetaR{i}' * v);
            
            % Ensure normalization of the perturbed stationary distribution
            dv = dv / sum(v) - v_norm * (sum(dv) / sum(v));
            
            % Reshape the derivative vector and map it to detection probabilities
            dV = reshape(dv, Ns, []);
            dM = E' * dV;
            d_thetaR{i} = foldVector(dM);
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

    % Reshape the result into a row vector
    Ps = reshape(M', 1, []);
end