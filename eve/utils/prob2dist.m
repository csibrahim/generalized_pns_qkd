function [d, dp, da] = prob2dist(p, alpha)
    %prob2dist: This function computes the distance (d) over which the signal 
    %            is attenuated, given the probability (p) and attenuation coefficient (alpha).
    %            It can also compute the partial derivatives of d with respect to p and alpha 
    %            if requested.
    %
    %            The distance is computed as d = -10 * log10(p) / alpha, which models the 
    %            inverse of the attenuation process.
    %
    % Inputs:
    %     p      - Probability of signal detection after attenuation
    %     alpha  - Attenuation coefficient
    %
    % Outputs:
    %     d      - Distance over which the signal is attenuated to the probability p
    %     dp     - Partial derivative of d with respect to p
    %     da     - Partial derivative of d with respect to alpha
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Compute the distance based on the probability and attenuation coefficient
    d = -10 * log10(p) ./ alpha;

    % If derivatives are requested, compute the partial derivatives
    if nargout > 1
        % Partial derivative with respect to probability (p)
        dp = d ./ (p .* log10(p) * log(10));

        % Partial derivative with respect to attenuation coefficient (alpha)
        da = -d ./ alpha;
    end

end