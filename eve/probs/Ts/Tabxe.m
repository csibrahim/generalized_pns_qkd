function [Ts, d_thetaE] = Tabxe(thetaA, thetaB, thetaE, a, b, x, e)
    %Tabxe: Computes the transition matrix and its derivatives for all
    %       intensities (lambdas), w.r.t. the parameters in thetaE. 
    %       The detection probabilities depend on Alice's and Bob's basis 
    %       choices (a, b), Alice's bit choice (x), and Eve's interception 
    %       flag (e)
    %
    % Inputs:
    %   thetaA - Alice’s parameter set [lambdas, alpha, dAB]
    %   thetaB - Bob’s parameter set [pa0, pa1, pc0, pc1, pd0, pd1, pe]
    %   thetaE - Eve’s parameter set [dAE, pEB, k, Delta]
    %   a      - Binary input, Alice's basis choice (0 or 1)
    %   b      - Binary input, Bob's basis choice (0 or 1)
    %   x      - Binary input, Alice's bit choice (0 or 1)
    %   e      - Binary input, Eve's interception flag 
    %            (1 if interception occurs, 0 otherwise)
    %
    % Outputs:
    %   Ts       - Combined transition matrix for all intensities (lambdas)
    %   d_thetaE - Cell array of partial derivatives of the transition matrix 
    %              with respect to the parameters in thetaE
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Retrieve lambdas and after-pulse probabilities from thetaA and thetaB
    lambdas = thetaA{1};  % Pulse intensities
    pa0 = thetaB{1};      % After-pulse probability for detector D0
    pa1 = thetaB{2};      % After-pulse probability for detector D1
    pa = [pa0, pa1];      % Group after-pulse probabilities

    Nl = length(lambdas); % Number of intensities

    if nargout > 1
        % Compute independent detection probabilities (without after-pulsing)
        [Ps, d_thetaE] = Pabxe(thetaA, thetaB, thetaE, a, b, x, e);

        % Preallocate for partial derivatives w.r.t. each intensity and
        % parameter in thetaE
        d_thetas = cell(Nl, numel(thetaE));
    else
        % If no derivatives are requested, just compute the independent 
        % detection probabilities (without after-pulsing)
        Ps = Pabxe(thetaA, thetaB, thetaE, a, b, x, e);
    end

    % Preallocate transition matrices for each intensity
    Ts = cell(Nl, 1);

    for i = 1:Nl
        
        % Select detection probabilities for the current intensity
        idx = false(1, length(Ps));
        idx((i - 1) + 1:Nl:length(idx)) = true;
        
        % Scale probabilities for the current intensity to ensure normalization
        Ps_i = Ps(idx) * Nl;

        if nargout > 1
            
            % Compute the transition matrix and its derivatives for lambda_i
            [Ts{i}, d_p00, d_p01, d_p10, d_p11] = Tabxel(Ps_i, pa);
            
            % Compute partial derivatives for each parameter in thetaE
            for j = 1:numel(d_thetaE)
                d_theta_j = d_thetaE{j}(idx) * Nl;
                d_thetas{i, j} = d_p00 * d_theta_j(1) + ...
                                 d_p01 * d_theta_j(2) + ...
                                 d_p10 * d_theta_j(3) + ...
                                 d_p11 * d_theta_j(4);
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
        for i = 1:numel(d_thetaE)
            d_thetaE{i} = horzcat(d_thetas{:, i}) / Nl;
            d_thetaE{i} = repmat(d_thetaE{i}, Nl, 1);
        end
    end
end











