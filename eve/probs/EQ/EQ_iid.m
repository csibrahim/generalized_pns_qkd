function [EQs, Qs] = EQ_iid(thetaA, thetaB, thetaE, both)
    %EQ_iid: Computes signal (Q) and error (EQ) probabilities for the i.i.d. 
    %        case across different pulse intensities. The function handles 
    %        both matching and non-matching basis configurations and supports 
    %        XOR detection or double-clicks based on the 'both' flag.
    %
    % Inputs:
    %   thetaA - Alice’s parameter set [lambdas, alpha, dAB]
    %   thetaB - Bob’s parameter set [pa0, pa1, pc0, pc1, pd0, pd1, pe]
    %   thetaE - Eve’s parameter set [dAE, pEB, k, Delta]
    %   both   - Logical flag (default: false):
    %            - false: Only XORs are considered signal events
    %            - true:  Includes ANDs as signal and error events
    %
    % Outputs:
    %   EQs - Error probabilities for each intensity and configuration
    %    Qs - Signal probabilities for each intensity and configuration
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    if nargin < 4
        both = false;  % Default to XOR detection
    end

    % Compute probabilities for non-matching basis configurations (01 and 10)
    [EQs01, Qs01] = EQ_iid_ab(thetaA, thetaB, thetaE, 0, 1, both);
    [EQs10, Qs10] = EQ_iid_ab(thetaA, thetaB, thetaE, 1, 0, both);

    % Average probabilities for non-matching bases
    EQs0 = (EQs01 + EQs10) / 2;
    Qs0 = (Qs01 + Qs10) / 2;

    % Compute probabilities for the matching basis configuration (11)
    [EQs1, Qs1] = EQ_iid_ab(thetaA, thetaB, thetaE, 1, 1, both);

    % Combine matching and non-matching probabilities
    EQs = [EQs0, EQs1] / 2;
    Qs = [Qs0, Qs1] / 2;
end

function [EQs, Qs] = EQ_iid_ab(thetaA, thetaB, thetaE, a, b, both)
    %EQ_iid_ab: Computes signal (Q) and error (EQ) probabilities for a 
    %           specific basis configuration (a, b) by marginalizing over 
    %           Alice's bit choice (x).
    
    % Compute probabilities for both x = 0 and x = 1
    [EQs0, Qs0] = EQ_iid_abx(thetaA, thetaB, thetaE, a, b, 0, both);
    [EQs1, Qs1] = EQ_iid_abx(thetaA, thetaB, thetaE, a, b, 1, both);

    % Average probabilities over Alice's bit choice (x)
    EQs = (EQs0 + EQs1) / 2;
    Qs = (Qs0 + Qs1) / 2;
end

function [EQs, Qs] = EQ_iid_abx(thetaA, thetaB, thetaE, a, b, x, both)
    %EQ_iid_abx: Computes signal (Q) and error (EQ) probabilities for a 
    %            specific basis configuration (a, b) and Alice's bit 
    %            choice (x).

    % Get detection probabilities for the given configuration
    Ps = Pabx(thetaA, thetaB, thetaE, a, b, x);

    Nl = size(Ps, 2) / 4;  % Number of intensities

    p0 = Ps(:, (Nl+1):2*Nl);   % prrobability of click at D0 only
    p1 = Ps(:, (2*Nl+1):3*Nl); % prrobability of click at D1 only

    % Define correct and error probabilities based on x
    if x == 0
        p_correct = p0;
        p_error   = p1;
    else
        p_correct = p1;
        p_error   = p0;
    end

    % Compute the total signal (correct + error)
    p_signal = p_correct + p_error;

    % If 'both' is true, include double-click to signal and error
    if both
        
        % Get double-click probability
        p_both = Ps(:, (3*Nl+1):4*Nl);

        % Add to both signal and error
        p_error  = p_error + p_both;
        p_signal = p_signal + p_both;
    end

    % Return the final signal and error probabilities
    EQs = p_error;
    Qs = p_signal;
end
