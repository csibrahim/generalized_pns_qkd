function [deltaQs, Qs] = deltaQ_iid(thetaR, thetaF, varR, varF, both)
    %deltaQ_iid: This function computes the probabilities of detection (Qs) and 
    %   the error probabilities (deltaQs) for the independent and identically 
    %   distributed (iid) case. The probabilities are computed based on whether 
    %   a click occurs at either of the detectors or from XOR (exclusive or) logic 
    %   between the detectors, depending on the 'both' flag.
    %
    %   If 'both' is true, the function computes Qs as the sum of P01 + P10 + P11, 
    %   which represents the detection event occurring at either or both detectors. 
    %   If 'both' is false, Qs is computed as the sum of P01 + P10, representing 
    %   the XOR of the detection events at the two detectors.
    %
    %   The function returns probabilities for both matching and non-matching 
    %   basis configurations, and across different pulse intensities.
    %
    % Inputs:
    %   thetaR - Array of random system parameters
    %   thetaF - Array of fixed system parameters
    %   varR   - Cell array of variable names corresponding to the random parameters
    %   varF   - Cell array of variable names corresponding to the fixed parameters
    %   both   - Logical flag. When true, the function considers clicks from both 
    %            detectors. When false, it computes the XOR detection probability 
    %            (default: false)
    %
    % Outputs:
    %   deltaQs - Error probabilities for each intensity and configuration
    %   Qs      - Detection probabilities for each intensity and configuration
    
    if nargin < 5
        both = false;  % Default: consider XOR detection
    end

    % Compute probabilities for the 01 and 10 basis configurations
    [deltaQs01, Qs01] = deltaQ_iid_ab(thetaR, thetaF, varR, varF, 0, 1, both);
    [deltaQs10, Qs10] = deltaQ_iid_ab(thetaR, thetaF, varR, varF, 1, 0, both);
    
    % Average the results for non-matching basis configurations (01 and 10)
    deltaQs0 = (deltaQs01 + deltaQs10) / 2;
    Qs0 = (Qs01 + Qs10) / 2;
    
    % Compute probabilities for the 11 (matching) basis configuration
    [deltaQs1, Qs1] = deltaQ_iid_ab(thetaR, thetaF, varR, varF, 1, 1, both);
    
    % Combine matching and non-matching probabilities
    deltaQs = [deltaQs0, deltaQs1] / 2;
    Qs = [Qs0, Qs1] / 2;
end

function [deltaQs, Qs] = deltaQ_iid_ab(thetaR, thetaF, varR, varF, a, b, both)
    %deltaQ_iid_ab: This function computes the probabilities for a given 
    %   basis configuration (a, b) in the iid case, by marginalizing over Alice's 
    %   bit choice.

    % Compute probabilities for both x = 0 and x = 1
    [deltaQs0, Qs0] = deltaQ_iid_abx(thetaR, thetaF, varR, varF, a, b, 0, both);
    [deltaQs1, Qs1] = deltaQ_iid_abx(thetaR, thetaF, varR, varF, a, b, 1, both);

    % Average the results over Alice's bit choice (x)
    deltaQs = (deltaQs0 + deltaQs1) / 2;
    Qs = (Qs0 + Qs1) / 2;
end

function [deltaQs, Qs] = deltaQ_iid_abx(thetaR, thetaF, varR, varF, a, b, x, both)
    %deltaQ_iid_abx: This function computes the error and detection probabilities
    %   for a specific basis configuration (a, b) and bit choice (x) in the iid case. 
    %   It uses the detection probabilities P from the Pabx function. Based on the 
    %   'both' flag, it computes either the XOR detection (P01 + P10) or includes 
    %   both detectors (P01 + P10 + P11).

    % Get detection probabilities for the given configuration
    Ps = Pabx(thetaR, thetaF, varR, varF, a, b, x);
    
    Nl = size(Ps, 2) / 4;  % Number of intensities

    % For x = 0, the correct and error terms are swapped
    if x == 0
        p_correct = Ps(:, (Nl+1):2*Nl);
        p_error = Ps(:, (2*Nl+1):3*Nl);
    else
        p_error = Ps(:, (Nl+1):2*Nl);
        p_correct = Ps(:, (2*Nl+1):3*Nl);
    end
    
    % Compute the total signal (correct + error)
    p_signal = p_correct + p_error;

    % If 'both' is true, include the case where both detectors click (P11)
    if both
        p_both = Ps(:, (3*Nl+1):4*Nl);
        p_error = p_error + p_both;
        p_signal = p_signal + p_both;
    end

    % Return the error probabilities and the total detection probabilities
    deltaQs = p_error;
    Qs = p_signal;
end
