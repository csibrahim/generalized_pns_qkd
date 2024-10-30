function varargout = Tabxe(thetaA, thetaB, thetaE, a, b, x, e)
    %Tabxe: This function computes the transition matrices for Bob's detectors 
    %       for all possible photon intensities (lambdas). It returns the combined 
    %       transition matrix and, if requested, the partial derivatives of the 
    %       transition matrix with respect to the parameters in thetaE.
    %       The detection probabilities are affected by Alice's and Bob's basis 
    %       choices (a, b), Alice's bit choice (x), and whether Eve intercepted 
    %       the pulse (e).
    %
    % Inputs:
    %   thetaA - Cell array of system parameters specific to the channel between 
    %            Alice and Bob, containing [lambdas, alpha, dAB].
    %   thetaB - Cell array of system parameters for Bob's detectors, containing 
    %            [pa0, pa1, pc0, pc1, pd0, pd1, pe].
    %   thetaE - Cell array of system parameters related to Eve’s eavesdropping, 
    %            containing [dAE, pEB, k, Delta].
    %   a      - Binary input, Alice's basis choice for the pulse (0 or 1).
    %   b      - Binary input, Bob's basis choice for the pulse (0 or 1).
    %   x      - Binary input, Alice's bit choice for the pulse (0 or 1).
    %   e      - Binary input, Eve's interception flag 
    %            (e=1 if she intercepted the pulse, e=0 otherwise).
    %
    % Outputs:
    %   varargout{1} - Combined transition matrix for all possible intensities (lambdas).
    %   varargout{2} - Cell array of partial derivatives of the transition matrix 
    %                  with respect to the parameters in thetaE.
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Initialize output cells
    varargout = cell(1, max(nargout, 1));

    % Retrieve lambdas and after-pulse probabilities from thetaA and thetaB
    lambdas = thetaA{1};
    pa0 = thetaB{1};
    pa1 = thetaB{2};
    pa = [pa0, pa1];

    Nl = length(lambdas);  % Number of intensities

    % Compute detection probabilities without after-pulsing
    [varargout{:}] = Pabxe(thetaA, thetaB, thetaE, a, b, x, e);
    Ps = varargout{1};

    Ts = cell(Nl, 1);  % Preallocate for transition matrices
    d_thetas = cell(Nl, numel(thetaE));  % Preallocate for derivatives

    for i = 1:Nl
        % Select detection probabilities for the current intensity
        idx = false(1, length(Ps));
        idx((i - 1) + 1:Nl:length(idx)) = true;
        
        % Scale probabilities for the current intensity to normalize to 1
        Ps_i = Ps(idx) * Nl;

        if nargout > 1
            % Retrieve derivatives of detection probabilities
            d_thetaE = varargout{2};

            % Compute the transition matrix and its derivatives
            [Ts{i}, d_p00, d_p01, d_p10, d_p11] = Tabxel(Ps_i, pa);
            
            % Calculate the partial derivatives for each parameter in thetaE
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
    varargout{1} = Ts;

    if nargout > 1
        % Concatenate and normalize the derivatives across intensities
        for i = 1:numel(d_thetaE)
            d_thetaE{i} = horzcat(d_thetas{:, i}) / Nl;
            d_thetaE{i} = repmat(d_thetaE{i}, Nl, 1);
        end
        varargout{2} = d_thetaE;
    end
end