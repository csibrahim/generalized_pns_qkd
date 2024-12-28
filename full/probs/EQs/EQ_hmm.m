function [EQs, Qs] = EQ_hmm(thetaR, thetaF, varR, varF, both)
    %EQ_hmm: Computes signal (Q) and error (EQ) probabilities for the HMM 
    %        case across different pulse intensities. The function handles 
    %        both matching and non-matching basis configurations and supports 
    %        XOR detection or double-clicks based on the 'both' flag.
    %
    % Inputs:
    %   thetaR  - Values of parameters listed in varR 
    %   thetaF  - Values of parameters listed in varF 
    %   varR    - Cell array of random variable names
    %   varF    - Cell array of fixed variable names
    %   both    - Logical flag (default: false):
    %             - false: Only XORs are considered signal events
    %             - true:  Includes ANDs as signal and error events
    %
    % Outputs:
    %   EQs - Error probabilities for each intensity and configuration
    %    Qs - Signal probabilities for each intensity and configuration
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    if nargin < 5
        both = 0;  % Default: consider XOR detection
    end

    % The emission matrix
    E = [1 0 0 0;  % p_00
         0 1 0 0;  % p_01
         0 0 1 0;  % p_10
         0 0 0 1;  % p_11
         0 1 0 0;  % p_0a
         0 0 1 0;  % p_a0
         0 0 0 1;  % p_1a
         0 0 0 1;  % p_a1
         0 0 0 1]; % p_aa

    % Number of states
    Ns = size(E, 1);

    % Compute the transition matrix
    Ts = T_full(thetaR, thetaF, varR, varF);

    % Compute the stationary distribution
    [v, ~] = eigs(Ts', 1, "largestreal");
    v = abs(v);           % Ensure non-negativity
    v_norm = v / sum(v);  % Normalize stationary distribution
    
    % Reshape the stationary vector to match the state space
    V = reshape(v_norm, Ns, []);

    % Map the stationary distribution to detection probabilities
    Ps = E' * V;

    % Number of intensities
    Nl = size(Ps, 2) / (3 * 2^2);

    
    Ps0 = Ps(:, 1:Nl*8);       % non-matching bases (a = 1, b = 0) and (a = 0, b = 1)
    Ps1 = Ps(:, (Nl*8+1):end); % matching bases (a = 1, b = 1)

    % Merge (a = 1, b = 0) with (a = 0, b = 1)
    Ps0 = Ps0(:, 1:4*Nl) + Ps0(:, 4*Nl+1:end);
    
    % Extract signal and error probabilities for matching and non-matching
    [EQs0, Qs0] = extract_EQs(Ps0, both);
    [EQs1, Qs1] = extract_EQs(Ps1, both);

    % Combine results for matching and non-matching bases
    EQs = [EQs0, EQs1];
    Qs = [Qs0, Qs1];
end

function [EQs, Qs] = extract_EQs(Ps, both)
    %extract_EQs: Computes signal (Q) and error (EQ) probabilities from a 
    %             given probability matrix. Handles XOR or simultaneous clicks 
    %             based on the 'both' flag.
    %
    % Inputs:
    %   Ps    - Probability matrix for the current basis configuration
    %   both  - Logical flag (default: false):
    %           - false: Only XORs are considered signal events
    %           - true:  Includes ANDs as signal and error events
    %
    % Outputs:
    %   EQs - Error probabilities for each intensity and configuration
    %    Qs - Signal probabilities for each intensity and configuration

    % Number of intensities
    Nl = size(Ps, 2) / (2^2);

    % Separate probabilities for x = 0 and x = 1
    Ps0 = Ps(:, 1:2*Nl);     % x = 0
    Ps1 = Ps(:, 2*Nl+1:end); % x = 1

    % Combine probabilities
    Ps1 = [Ps1(1, :); Ps1(3, :); Ps1(2, :); Ps1(4, :)];

    % Marginalize over Eve's interception flag (e)
    Ps = Ps0 + Ps1;

    % Merge Alice and Eve's bits
    Ps = Ps(:, 1:Nl) + Ps(:, Nl+1:end);

    % Now Ps has four rows [no-click, correct, error, double-click] and Nl columns
     correct = Ps(2,:);
       error = Ps(3,:);
    double_click = Ps(4,:);

     Qs = correct+error;
    EQs = error;

    % If both is true, add double-clicks to signals and errors
    if both
        Qs = Qs + double_click;
        EQs = EQs + double_click;
    end
end