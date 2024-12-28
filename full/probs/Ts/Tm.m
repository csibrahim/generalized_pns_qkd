function [Ts, d_thetaR] = Tm(thetaR, thetaF, varR, varF, m)
    %Tm: Computes the transition matrix and its derivatives for matching
    %    or non-matching bases (m) w.r.t. the parameters in thetaR
    %    marginalzing over Alice's bit choise (x) and Eve's
    %    interception flag (e)
    %
    % Inputs:
    %   thetaR - Values of parameters listed in varR 
    %   thetaF - Values of parameters listed in varF 
    %   varR   - Cell array of random variable names
    %   varF   - Cell array of fixed variable names
    %   m      - Basis configuration: 0 = non-matching, 1 = matching
    %
    % Outputs:
    %   Ts       - Transition matrix for matching or non-matching bases (m)
    %              marginalzing over Alice's bit choise (x) and Eve's 
    %              interception flag (e)
    %   d_thetaR - Cell array of partial derivatives of the transition matrix 
    %              with respect to the parameters in thetaR
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    a = 1; % Assume Alice always prepares in this basis
    b = m; % Assume Bob always measures in this basis

    % Concatenate random and fixed variables for easier access
    vars = [varR(:)' varF(:)'];
    theta = [thetaR(:)' thetaF(:)'];

    % Identify the indices of key parameters in the concatenated list
    [~, ind_lambdas] = ismember('lambdas', vars);
    [~, ind_pa0] = ismember('pa0', vars);
    [~, ind_pa1] = ismember('pa1', vars);

    % Retrieve lambdas and after-pulse probabilities
    lambdas = theta{ind_lambdas}; % Pulse intensities
    pa0 = theta{ind_pa0};         % After-pulse probability for detector D0
    pa1 = theta{ind_pa1};         % After-pulse probability for detector D1
    pa = [pa0, pa1];              % Group after-pulse probabilities

    Nl = length(lambdas);  % Number of intensities

    if nargout > 1
        % Compute independent detection probabilities (without after-pulsing)
        [Ps, d_thetaR] = Pab(thetaR, thetaF, varR, varF, a, b);
        
        % Preallocate storage for derivatives
        d_thetas = cell(Nl, numel(thetaR));

        % Recompute indices for lambdas, pa0, and pa1 in the random variables
        [~, ind_lambdas] = ismember('lambdas', varR);
        [~, ind_pa0] = ismember('pa0', varR);
        [~, ind_pa1] = ismember('pa1', varR);

        % Adjust indices for pa0 and pa1 if lambdas come earlier
        if ind_lambdas && ind_pa0 > ind_lambdas
            ind_pa0 = ind_pa0 + Nl - 1;
        end
        if ind_lambdas && ind_pa1 > ind_lambdas
            ind_pa1 = ind_pa1 + Nl - 1;
        end
    else
        % If no derivatives are requested, just compute the independent 
        % detection probabilities (without after-pulsing)
        Ps = Pab(thetaR, thetaF, varR, varF, a, m);
    end

    % Preallocate transition matrices for each intensity
    Ts = cell(Nl, 1);

    for i = 1:Nl

        % Select detection probabilities for the current intensity
        idx = false(1, length(Ps));
        idx((i - 1) + 1:Nl:length(idx)) = true;
        
        % Scale probabilities for the current intensity to ensure normalization.
        Ps_i = Ps(idx) * Nl;

        if nargout > 1
            
            % Compute the transition matrix and its derivatives for lambda_i
            [Ts{i}, d_p00, d_p01, d_p10, d_p11, d_pa0, d_pa1] = Tabxel(Ps_i, pa);

            % Compute partial derivatives for each parameter in thetaR
            for j = 1:numel(d_thetaR)
                if j == ind_pa0
                    d_thetas{i, j} = d_pa0;
                elseif j == ind_pa1
                    d_thetas{i, j} = d_pa1;
                else
                    d_theta_j = d_thetaR{j}(idx) * Nl;
                    d_thetas{i, j} = d_p00 * d_theta_j(1) + d_p01 * d_theta_j(2) + d_p10 * d_theta_j(3) + d_p11 * d_theta_j(4);
                end
            end
        else
            % Compute only the transition matrix for the current intensity
            Ts{i} = Tabxel(Ps_i, pa);
        end
    end

    % Concatenate and normalize the transition matrices across intensities
    Ts = horzcat(Ts{:}) / Nl;
    Ts = repmat(Ts, Nl, 1);

    if nargout > 1
        % Concatenate and normalize the derivatives across intensities
        for i = 1:numel(d_thetaR)
            d_thetaR{i} = horzcat(d_thetas{:, i}) / Nl;
            d_thetaR{i} = repmat(d_thetaR{i}, Nl, 1);
        end
    end
end