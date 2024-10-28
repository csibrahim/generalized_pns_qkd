function [Ts, d_thetaR] = Tabxe(thetaR, thetaF, varR, varF, a, b, x, e)
    %Tabxe: This function computes the transition matrices for Bob's detectors 
    %       for all possible photon intensities (lambdas). It also returns the 
    %       derivatives of the transition matrices with respect to the random variables
    %       specified in varR, if requested. The detection probabilities are affected
    %       by Alice's and Bob's basis choices (a, b), Alice's bit choice (x), and 
    %       whether Eve intercepted the pulse (e).

    % Inputs:
    %   thetaR - Array of system parameters treated as random variables 
    %            (e.g., lambdas, pa0, pa1, etc.)
    %   thetaF - Array of system parameters treated as fixed variables
    %   varR   - Cell array of variable names for the random system parameters
    %   varF   - Cell array of variable names for the fixed system parameters
    %   a      - Binary input, Alice's basis choice for the pulse (0 or 1)
    %   b      - Binary input, Bob's basis choice for the pulse (0 or 1)
    %   x      - Binary input, Alice's bit choice for the pulse (0 or 1)
    %   e      - Binary input, Eve's interception flag 
    %            (e=1 -> she intercepted this pulse, e=0 -> she did not intercept it)
    %
    % Outputs:
    %   Ts      - Combined transition matrix for all possible intensities (lambdas)
    %   d_thetaR - (Optional) Cell array of partial derivatives of the transition 
    %              matrix with respect to the random variables in varR
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Concatenate random and fixed variables for easier access
    vars = [varR(:)' varF(:)'];
    theta = [thetaR(:)' thetaF(:)'];

    % Identify the indices of key parameters in the concatenated list
    [~, ind_lambdas] = ismember('lambdas', vars);
    [~, ind_pa0] = ismember('pa0', vars);
    [~, ind_pa1] = ismember('pa1', vars);

    % Retrieve the values for lambdas and after-pulse probabilities
    lambdas = theta{ind_lambdas};
    pa = [theta{ind_pa0}, theta{ind_pa1}];

    Nl = length(lambdas);  % Number of intensities

    if nargout > 1
        % Compute detection probabilities and their partial derivatives
        [Ps, d_thetaR] = Pabxe(thetaR, thetaF, varR, varF, a, b, x, e);
        d_thetas = cell(Nl, numel(thetaR));  % Preallocate storage for derivatives

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
        % If no derivatives are requested, just compute the detection probabilities
        Ps = Pabxe(thetaR, thetaF, varR, varF, a, b, x, e);
    end

    Ts = cell(Nl, 1);  % Preallocate for transition matrices

    for i = 1:Nl
        % Select detection probabilities for the current intensity
        idx = false(1, length(Ps));
        idx((i - 1) + 1:Nl:length(idx)) = true;
        
        % Scale probabilities for the current intensity
        Ps_i = Ps(idx) * Nl;

        if nargout > 1
            % Compute the transition matrix and its derivatives
            [Ts{i}, d_p00, d_p01, d_p10, d_p11, d_pa0, d_pa1] = Tabxel(Ps_i, pa);

            % Compute the derivatives for each random parameter
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