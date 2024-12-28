function varargout = P(thetaR, thetaF, varR, varF)
    %P: Computes detection probabilities and thier derivatives w.r.t. the 
    %   parameters in thetaR using the appropriate model based on the system
    %   parameters. If either after-pulse probabilities (pa0, pa1) are 
    %   greater than zero, the hidden Markov model (HMM) based function
    %   P_hmm is called. Otherwise, the independent and identically 
    %   distributed (i.i.d.) function P_iid is used.
    %
    % Inputs:
    %   thetaR - Values of parameters listed in varR 
    %   thetaF - Values of parameters listed in varF 
    %   varR   - Cell array of random variable names
    %   varF   - Cell array of fixed variable names
    %
    % Outputs:
    %   Ps       - Combined detection probabilities for all intensities 
    %              (lambdas), matching and non-matching bases marginalzed 
    %              over Alice's bit choise (x) and Eve's interception flag (e)
    %   d_thetaR - Cell array of partial derivatives of detection probabilities 
    %              with respect to the parameters in thetaR
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Initialize output cells
    varargout = cell(1, max(1, nargout));

    % Concatenate random and fixed variables for easy access
    vars = [varR(:)' varF(:)'];
    theta = [thetaR(:)' thetaF(:)'];

    % Identify the indices for the after-pulse probabilities pa0 and pa1
    [~, ind_pa] = ismember({'pa0', 'pa1'}, vars);
    
    % Retrieve the values of pa0 and pa1
    pa = [theta{ind_pa}];

    % Select the appropriate model based on after-pulse probabilities
    if any(pa > 0)
        % Call the HMM-based function if after-pulse probabilities are non-zero
        [varargout{:}] = P_hmm(thetaR, thetaF, varR, varF);
    else
        % Call the iid-based function if there are no after-pulses
        [varargout{:}] = P_iid(thetaR, thetaF, varR, varF);
    end
end
