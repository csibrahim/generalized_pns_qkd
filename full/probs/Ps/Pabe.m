function [Ps, d_thetaR] = Pabe(thetaR, thetaF, varR, varF, a, b, e)
    %Pabe: Computes detection probabilities and their derivatives for all
    %      intensities (lambdas) and marginalized over Alice's bit choice 
    %      (x) w.r.t. the parameters in thetaR. The detection probabilities 
    %      depend on Alice's and Bob's basis choices (a, b), and Eve's 
    %      interception flag (e)
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
    %              (lambdas) marginalized over Alice's bit choice (x)
    %   d_thetaR - Cell array of partial derivatives of the detection 
    %              probabilities with respect to the parameters in thetaR
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Check if partial derivatives are requested
    if nargout > 1

        % Compute detection probabilities and derivatives for x = 0 and x = 1
        [Ps0, d_thetaR0] = Pabxe(thetaR, thetaF, varR, varF, a, b, 0, e);
        [Ps1, d_thetaR1] = Pabxe(thetaR, thetaF, varR, varF, a, b, 1, e);
        
        % Preallocate space for derivatives
        d_thetaR = cell(1, numel(d_thetaR0));
    
        % Compute averaged partial derivatives for each parameter in thetaR
        for i = 1:numel(d_thetaR)
            d_thetaR{i} = (d_thetaR0{i} + d_thetaR1{i}) / 2;
        end
        
    else

        % Compute detection probabilities for x = 0 and x = 1
        Ps0 = Pabxe(thetaR, thetaF, varR, varF, a, b, 0, e);
        Ps1 = Pabxe(thetaR, thetaF, varR, varF, a, b, 1, e);

    end

    % Marginalize over Alice's bit choice (x)
    Ps = (Ps0 + Ps1) / 2;

end