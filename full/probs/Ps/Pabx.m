function [Ps, d_thetaR] = Pabx(thetaR, thetaF, varR, varF, a, b, x)
    %Pabx: Computes the detection probabilities and their derivatives for all
    %      intensities (lambdas) and Eve's interception flag (e) w.r.t. 
    %      the parameters in thetaR. The detection probabilities depend on
    %      Alice's and Bob's basis choices (a, b), Alice's bit choice (x)
    %
    % Inputs:
    %   thetaR - Values of parameters listed in varR 
    %   thetaF - Values of parameters listed in varF 
    %   varR   - Cell array of random variable names
    %   varF   - Cell array of fixed variable names
    %   a      - Binary input, Alice's basis choice (0 or 1)
    %   b      - Binary input, Bob's basis choice (0 or 1)
    %   x      - Binary input, Alice's bit choice (0 or 1)
    %
    % Outputs:
    %   Ps       - Combined detection probabilities for all intensities 
    %              (lambdas) and Eve's interception flag (e)
    %   d_thetaR - Cell array of partial derivatives of the detection 
    %              probabilities with respect to the parameters in thetaR
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Combine random and fixed variables
    vars = [varR(:)' varF(:)'];
    theta = [thetaR(:)' thetaF(:)'];

    % Find the index of Delta in the parameter list
    [~, ind_Delta] = ismember('Delta', vars);
    Delta = theta{ind_Delta};

    % If derivatives are requested
    if nargout > 1
        % Compute detection probabilities and derivatives for e = 0 and e = 1
        [Ps0, d_thetaR0] = Pabxe(thetaR, thetaF, varR, varF, a, b, x, 0);  % No interception
        [Ps1, d_thetaR1] = Pabxe(thetaR, thetaF, varR, varF, a, b, x, 1);  % Interception

        % Preallocate space for derivatives
        d_thetaR = cell(1, numel(d_thetaR0));
        
        % Combine derivatives for each parameter weighted by Delta
        for i = 1:numel(d_thetaR)
            d_thetaR{i} = Delta .* d_thetaR1{i} + (1 - Delta) .* d_thetaR0{i};
        end

        % If Delta is one of the random variables, compute its derivative
        [~, ind_Delta_in_varR] = ismember('Delta', varR);
        
        if ind_Delta_in_varR
            [in_varR, ind_lambda] = ismember('lambdas', varR);

            % Adjust Delta's index if lambdas come earlier in the variable list
            if in_varR && ind_lambda < ind_Delta_in_varR
                lambdas = thetaR{ind_lambda};
                Nl = size(lambdas, 2);
                ind_Delta_in_varR = ind_Delta_in_varR + Nl - 1;
            end

            % Derivative w.r.t. Delta
            d_thetaR{ind_Delta_in_varR} = Ps1 - Ps0;
        end

    else
        % Compute only the detection probabilities for e = 0 and e = 1
        Ps0 = Pabxe(thetaR, thetaF, varR, varF, a, b, x, 0);
        Ps1 = Pabxe(thetaR, thetaF, varR, varF, a, b, x, 1);
    end

    % Marginalize e using Delta
    Ps = Delta .* Ps1 + (1 - Delta) .* Ps0;

end