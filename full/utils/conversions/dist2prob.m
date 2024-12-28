function [p, dd, da] = dist2prob(d, alpha)
    %dist2prob: Computes the probability (p) from distance (d) and 
    %           attenuation coefficient (alpha). Optionally computes 
    %           partial derivatives with respect to d (dd) and alpha (da).
    %
    % Inputs:
    %     d      - Distance over which the signal is attenuated
    %     alpha  - Attenuation coefficient
    %
    % Outputs:
    %     p      - Probability of signal detection after distance (d) attenuation (alpha)
    %     dd     - Partial derivative of p w.r.t. distance (d)
    %     da     - Partial derivative of p w.r.t. attenuation (alpha)
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Compute the probability p based on the distance and attenuation
    p = 10.^(-alpha .* d / 10);

    % Compute partial derivatives if requested
    if nargout > 1
        % Common factor for the derivatives
        dp_de = -log(10) * p / 10;

        % Partial derivative w.r.t. distance (d)
        dd = alpha .* dp_de;

        % Partial derivative w.r.t. attenuation coefficient (alpha)
        da = d .* dp_de;
    end

end
