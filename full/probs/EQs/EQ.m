function [EQs, Qs, Deltas] = EQ(thetaR, thetaF, varR, varF, both)
    %EQ: Computes signal (Q) and error (EQ) probabilities probabilities
    %    based on whether the system follows an independent and identically 
    %    distributed (i.i.d.) model or a hidden Markov model (HMM). If 
    %    either after-pulse probabilities (pa0, pa1) are greater than zero,
    %    the hidden Markov model (HMM) based function EQ_hmm is called. 
    %    Otherwise, the independent and identically distributed (i.i.d.) 
    %    function EQ_iid is used.
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
    %    Qs - Detection probabilities for each intensity and configuration
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    if nargin < 5
        both = false;  % Default: consider XOR detection
    end

    if(~iscell(thetaR))
        thetaR = theta2cell(thetaR, varR);
    end

    % Concatenate random and fixed variables into a single array
    vars = [varR(:)' varF(:)'];
    thetas = [thetaR(:)' thetaF(:)'];

    % Extract number of intensities (Nl)
    [~, ind_lambdas] = ismember('lambdas', vars);
    Nl = size(thetas{ind_lambdas}, 2);

    % Extract Eve's interception rate Delta
    [~, ind_Delta] = ismember('Delta', vars);
    Deltas = thetas{ind_Delta}; 

    % Extract after-pulse probabilities (pa0, pa1)
    [~, ind_pa] = ismember({'pa0', 'pa1'}, vars);
    pa = [thetas{ind_pa}];

    % Determine if the HMM model should be used (non-zero after-pulse)
    if any(pa > 0) 

        % HMM case

        % Concatenate the random parameters into a single array
        thetaR = [thetaR{:}];

        n = size(thetaR,1);  % Number of samples

        % Initialize outputs for error and signal probabilities
        EQs = zeros(n, 2 * Nl);
        Qs = zeros(n, 2 * Nl);

        % Loop through each sample
        for i = 1:n
            % Convert thetaR(i, :) into cell format for the HMM case
            thetaR_i = theta2cell(thetaR(i, :), varR);

            % Compute the HMM probabilities
            [EQs(i, :), Qs(i, :)] = EQ_hmm(thetaR_i, thetaF, varR, varF, both);

            % Display progress if there are multiple samples
            if n > 2
                progress(i, n, 'Computing error and signal probabilities (Q and EQ)');
            end
        end
        
    else
        % iid case
        [EQs, Qs] = EQ_iid(thetaR, thetaF, varR, varF, both);
    end
end