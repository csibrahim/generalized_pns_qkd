function [thetaP, thetaR, thetaF, varR, varF] = processParams(thetas, true_thetas, varF, varR, sigma)
    %processParams: This function processes the input parameters, separating them into 
    %                fixed and random variables. For each random variable, it computes 
    %                the corresponding prior distribution parameters (alphas, betas, 
    %                upper and lower bounds) based on the distribution type (Gamma or Beta).
    %
    % Inputs:
    %     thetas      - Cell array containing the grouped variables: thetaA, thetaB, thetaE
    %     true_thetas - Cell array of true values for the parameters (ground truth)
    %     varF        - Cell array of fixed variable names
    %     varR        - Cell array of random variable names
    %     sigma       - Standard deviation for the prior distributions (default: 0)
    %
    % Outputs:
    %     thetaP      - Prior parameters for the random variables (alphas, betas, ub, lb)
    %     thetaR      - True values of the random variables
    %     thetaF      - Values of the fixed variables
    %     varR        - Updated random variable names
    %     varF        - Updated fixed variable names
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    if nargin < 5
        sigma = 0;  % Default sigma value if not provided
    end
    
    mus = struct();  % Structure to store the means of variables

    % Unpack theta parameters into respective groups (thetaA, thetaB, thetaE)
    [thetaA, thetaB, thetaE] = deal(thetas{:});
    [mus.lambdas, mus.alpha, mus.dAB] = deal(thetaA{:});
    [mus.pa0, mus.pa1, mus.pc0, mus.pc1, mus.pd0, mus.pd1, mus.pe] = deal(thetaB{:});
    [mus.dAE, mus.pEB, mus.k, mus.Delta] = deal(thetaE{:});

    % Get maximum value of 'k' for Eve's parameters
    k_max = get_k_max(thetaA, thetaB);

    % Initialize ground truth (gt) struct
    gt = struct();
    [thetaA, thetaB, thetaE] = deal(true_thetas{:});
    [gt.lambdas, gt.alpha, gt.dAB] = deal(thetaA{:});
    [gt.pa0, gt.pa1, gt.pc0, gt.pc1, gt.pd0, gt.pd1, gt.pe] = deal(thetaB{:});
    [gt.dAE, gt.pEB, gt.k, gt.Delta] = deal(thetaE{:});

    % Map variable names to distribution types
    distType = containers.Map();
    gammaVars = {'lambdas', 'alpha', 'dAB'};
    for i = 1:numel(gammaVars)
        distType(gammaVars{i}) = 'gamma';
    end
    betaVars = {'pa0', 'pa1', 'pc0', 'pc1', 'pd0', 'pd1', 'pe'};
    for i = 1:numel(betaVars)
        distType(betaVars{i}) = 'beta';
    end
    eveVars = {'dAE', 'pEB', 'k', 'Delta'};
    for i = 1:numel(eveVars)
        distType(eveVars{i}) = 'eve';
    end

    % Collect fixed parameters
    thetaF = cell(1, numel(varF));
    for i = 1:numel(varF)
        varName = varF{i};
        thetaF{i} = mus.(varName);
    end

    % Initialize prior parameter arrays
    thetaR = cell(1, numel(varR));  % True values of random variables
    alphas = cell(1, numel(varR));  % Alpha parameters for priors
    betas = cell(1, numel(varR));   % Beta parameters for priors
    ub = cell(1, numel(varR));      % Upper bounds for priors
    lb = cell(1, numel(varR));      % Lower bounds for priors

    % Process random variables
    for i = 1:numel(varR)
        varName = varR{i};
        thetaR{i} = gt.(varName);  % Assign ground truth values

        mu = mus.(varName);  % Mean value for the variable

        % Assign appropriate prior parameters based on distribution type
        switch distType(varName)
            case 'gamma'
                [alpha_param, beta_param, ub_param, lb_param] = gamma_parameters(mu, (mu * sigma).^2);
            case 'beta'
                [alpha_param, beta_param, ub_param, lb_param] = beta_parameters(mu, (mu * sigma).^2);
            case 'eve'
                switch varName
                    case 'dAE'
                        alpha_param = 1; beta_param = 2; ub_param = mus.dAB; lb_param = 0;
                    case 'pEB'
                        alpha_param = 1; beta_param = 1; ub_param = 1; lb_param = 0;
                    case 'k'
                        alpha_param = 1; beta_param = 2 / (k_max - 1); ub_param = inf; lb_param = 1;
                    case 'Delta'
                        alpha_param = 2; beta_param = 1; ub_param = 1; lb_param = 0;
                end
        end

        % Store prior parameters
        alphas{i} = alpha_param;
        betas{i} = beta_param;
        ub{i} = ub_param;
        lb{i} = lb_param;
    end

    % Pack prior parameters into a cell array
    thetaP = {alphas, betas, ub, lb};
end