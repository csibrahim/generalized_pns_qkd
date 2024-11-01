function [deltaQs, Qs, Deltas] = deltaQ(thetaR, thetaF, varR, varF, both)
    %deltaQ: This function computes the error probabilities (deltaQs) and 
    %        detection probabilities (Qs) based on whether the system follows an 
    %        independent and identically distributed (iid) model or a hidden Markov 
    %        model (HMM). The choice depends on the after-pulse probabilities (pa).
    %
    %        Additionally, the function extracts Deltas, which represent the 
    %        proportion of Eve's interceptions from the parameter set.
    %
    %        If the after-pulse probabilities (pa0, pa1) are greater than 0, the HMM 
    %        case is considered, and the function computes the probabilities based 
    %        on the hidden Markov model. Otherwise, the iid model is used.
    %
    %        The function computes these values for both matching and non-matching 
    %        basis configurations and across different pulse intensities.
    %
    % Inputs:
    %   thetaR - Array of random system parameters
    %   thetaF - Array of fixed system parameters
    %   varR   - Cell array of variable names corresponding to the random parameters
    %   varF   - Cell array of variable names corresponding to the fixed parameters
    %   both   - Logical flag, when true, the function considers clicks from both 
    %            detectors, when false it considers XOR (default: false)
    %
    % Outputs:
    %   deltaQs - Error probabilities for each intensity and configuration
    %   Qs      - Detection probabilities for each intensity and configuration
    %   Deltas  - Proportion of Eve's interceptions extracted from the parameters
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

    % Extract the number of intensities (Nl)
    [~, ind_lambdas] = ismember('lambdas', vars);
    Nl = size(thetas{ind_lambdas}, 2);

    % Extract the interception proportion (Delta), representing Eve's interception
    [~, ind_Delta] = ismember('Delta', vars);
    Deltas = thetas{ind_Delta}; 

    % Extract the after-pulse probabilities (pa0, pa1)
    [~, ind_pa] = ismember({'pa0', 'pa1'}, vars);
    pa = [thetas{ind_pa}];

    % Check if the system should follow the HMM model (non-zero after-pulse)
    if any(pa > 0)
        % HMM case

        % Concatenate the random parameters into a single array
        thetaR = [thetaR{:}];

        n = size(thetaR,1);  % Number of samples

        % Initialize output arrays for the error and detection probabilities
        deltaQs = zeros(n, 2 * Nl);
        Qs = zeros(n, 2 * Nl);

        % Loop through each sample
        for i = 1:n
            % Convert thetaR(i, :) into cell format for the HMM case
            thetaR_i = theta2cell(thetaR(i, :), varR);

            % Compute the HMM probabilities
            [deltaQs(i, :), Qs(i, :)] = deltaQ_hmm(thetaR_i, thetaF, varR, varF, both);

            % Update the progress indicator
            if n > 2
                progress(i, n, 'Computing QBERs (delta) and Gains (Q)');
            end
        end
        
    else
        % iid case
        [deltaQs, Qs] = deltaQ_iid(thetaR, thetaF, varR, varF, both);
    end
end