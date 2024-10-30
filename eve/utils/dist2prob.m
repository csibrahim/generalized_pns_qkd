function [p, dd, da] = dist2prob(d, alpha)
    %dist2prob: Computes the probability p given the distance (d) and 
    %            attenuation coefficient (alpha). The function can also compute the 
    %            partial derivatives of p with respect to the distance (dd) and 
    %            attenuation coefficient (da) if requested.
    %
    % Inputs:
    %     d      - Distance over which the signal is attenuated
    %     alpha  - Attenuation coefficient
    %
    % Outputs:
    %     p      - Probability of signal detection after attenuation
    %     dd     - Partial derivative of p with respect to distance (d)
    %     da     - Partial derivative of p with respect to the attenuation (alpha)
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Compute the probability p based on the distance and attenuation
    p = 10.^(-alpha .* d / 10);

    % If derivatives are requested, compute the partial derivatives
    if nargout > 1
        % Common factor in the derivatives
        dp_de = -log(10) * p / 10;

        % Partial derivative with respect to distance (d)
        dd = alpha .* dp_de;

        % Partial derivative with respect to attenuation coefficient (alpha)
        da = d .* dp_de;
    end

end