function [varR, varF] = processVars(thetas, varF, varR, sigma)
    %processVars: This function processes random and fixed variables by moving certain
    %              variables from the random set (varR) to the fixed set (varF) based
    %              on their prior means and distribution type (Gamma or Beta). Additionally,
    %              it adjusts the variables when the prior variance is zero.
    %
    % Inputs:
    %     thetas  - Cell array containing the grouped variables thetaA, thetaB, and thetaE
    %     varF    - Cell array of fixed variable names
    %     varR    - Cell array of random variable names
    %     sigma   - Prior variance; if zero, all non-Eve variables are moved to varF (default: 0)
    %
    % Outputs:
    %     varR    - Updated random variable names
    %     varF    - Updated fixed variable names
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    if (nargin < 4)
        sigma = 0;
    end
    
    mus = struct();  % Structure to hold the means of variables

    % Unpack the theta parameters into their respective groups
    [thetaA, thetaB, thetaE] = deal(thetas{:});
    
    % Assign thetaA, thetaB, and thetaE to their respective variables in `mus`
    [mus.lambdas, mus.alpha, mus.dAB] = deal(thetaA{:});
    [mus.pa0, mus.pa1, mus.pc0, mus.pc1, mus.pd0, mus.pd1, mus.pe] = deal(thetaB{:});
    [mus.dAE, mus.pEB, mus.k, mus.Delta] = deal(thetaE{:});

    % Remove Gamma-distributed variables with mean = 0 from varR to varF
    isGamma = ismember(varR, {'alpha', 'dAB'});  % Identify Gamma-distributed variables
    gammaVars = varR(isGamma);

    for i = 1:numel(gammaVars)
        name = gammaVars{i};
        [~, ind] = ismember(name, varR);  % Find the index in varR
        
        % Move the variable to varF if its mean is 0
        if (ind && mus.(name) == 0)
            varR(ind) = [];
            varF = [varF name];
        end
    end
    
    % Remove Beta-distributed variables with mean = 0 or 1 from varR to varF
    isBeta = ismember(varR, {'pa0', 'pa1', 'pc0', 'pc1', 'pd0', 'pd1', 'pe'});  % Identify Beta-distributed variables
    betaVars = varR(isBeta);

    for i = 1:numel(betaVars)
        name = betaVars{i};
        [~, ind] = ismember(name, varR);  % Find the index in varR

        % Move the variable to varF if its mean is 0 or 1
        if ind && (mus.(name) == 0 || mus.(name) == 1)
            varR(ind) = [];
            varF = [varF name];
        end
    end

    % If prior variance is set to 0, keep only Eve's parameters in varR
    if (sigma == 0)
        eveVars = {'dAE', 'k', 'pEB', 'Delta'};  % Eve's parameters
        isEveVar = ismember(varR, eveVars);

        % Move non-Eve variables from varR to varF
        varF = [varF varR(~isEveVar)];
        varR(~isEveVar) = [];
    end
    
end