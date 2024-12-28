function [Ps, d_thetaR] = Pab(thetaR, thetaF, varR, varF, a, b)
    %Pab: Computes the detection probabilities and their derivatives for all
    %     intensities (lambdas), marginalized over Eve's interception 
    %     flag (e), and Alice's bit choice (x) w.r.t. the parameters in 
    %     thetaR. The detection probabilities depend on Alice's and Bob's 
    %     basis choices (a, b)
    %
    % Inputs:
    %   thetaR - Values of parameters listed in varR 
    %   thetaF - Values of parameters listed in varF 
    %   varR   - Cell array of random variable names
    %   varF   - Cell array of fixed variable names
    %   a      - Binary input, Alice's basis choice (0 or 1)
    %   b      - Binary input, Bob's basis choice (0 or 1)
    %
    % Outputs:
    %   Ps       - Combined detection probabilities for all intensities 
    %              (lambdas), marginalized over Eve's interception flag (e) 
    %              and Alice's bit choice (x)
    %   d_thetaR - Cell array of partial derivatives of detection probabilities
    %              with respect to the parameters in thetaR
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Check if partial derivatives are requested
    if nargout > 1
        % Compute detection probabilities and derivatives for x = 0 and x = 1
        [Ps0, d_thetaR0] = Pabx(thetaR, thetaF, varR, varF, a, b, 0);
        [Ps1, d_thetaR1] = Pabx(thetaR, thetaF, varR, varF, a, b, 1);
        
        % Preallocate space for derivatives
        d_thetaR = cell(1, numel(d_thetaR0));
    
        % Average the derivatives for x = 0 and x = 1
        for i = 1:numel(d_thetaR)
            d_thetaR{i} = (d_thetaR0{i} + d_thetaR1{i}) / 2;
        end

    else
        % If no derivatives are requested, only compute probabilities
        Ps0 = Pabx(thetaR, thetaF, varR, varF, a, b, 0);
        Ps1 = Pabx(thetaR, thetaF, varR, varF, a, b, 1);
    end

    % Marginalize over Alice's bit choice (x)
    Ps = (Ps0 + Ps1) / 2;

end