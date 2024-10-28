function varargout = P(thetaR, thetaF, varR, varF)
    %P: This function computes the detection probabilities by selecting the 
    %   appropriate model based on the system parameters. If the after-pulse 
    %   probabilities (pa0, pa1) are greater than zero, it calls the hidden 
    %   Markov model (HMM) based function `P_hmm`. Otherwise, it assumes the 
    %   independent and identically distributed (iid) case and calls `P_iid`.
    %
    % Inputs:
    %   thetaR - Array of system parameters treated as random variables
    %   thetaF - Array of system parameters treated as fixed variables
    %   varR   - Cell array of variable names corresponding to the random parameters
    %   varF   - Cell array of variable names corresponding to the fixed parameters
    %
    % Outputs:
    %   varargout - Depending on the number of output arguments, this function
    %               returns detection probabilities or their partial derivatives 
    %               based on the selected model (iid or HMM)
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    varargout = cell(1, max(1, nargout));

    % Concatenate random and fixed variables for easy access
    vars = [varR(:)' varF(:)'];
    theta = [thetaR(:)' thetaF(:)'];

    % Identify the indices for the after-pulse probabilities pa0 and pa1
    [~, ind_pa] = ismember({'pa0', 'pa1'}, vars);
    pa = [theta{ind_pa}];  % Retrieve the values of pa0 and pa1

    % Select the appropriate model based on the after-pulse probabilities
    if any(pa > 0)
        % If after-pulse probabilities are greater than zero, use the HMM model
        [varargout{:}] = P_hmm(thetaR, thetaF, varR, varF);
    else
        % If no after-pulses, assume the iid case
        [varargout{:}] = P_iid(thetaR, thetaF, varR, varF);
    end
end
