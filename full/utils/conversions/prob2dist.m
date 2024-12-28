function [d, dp, da] = prob2dist(p, alpha)
    %prob2dist: Computes the distance (d) over which a signal attenuates
    %           to a given probability (p), using the attenuation
    %           coefficient (alpha). Optionally computes partial derivatives
    %           of d with respect to p and alpha.
    %
    % Inputs:
    %     p      - Probability of signal detection after attenuation
    %     alpha  - Attenuation coefficient
    %
    % Outputs:
    %     d      - Distance corresponding to the given probability (p) and attenuation (alpha)
    %     dp     - Partial derivative of d w.r.t. p
    %     da     - Partial derivative of d w.r.t. alpha
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Compute distance based on probability and attenuation coefficient
    d = -10 * log10(p) ./ alpha;

    % Compute partial derivatives if requested
    if nargout > 1
        % Partial derivative w.r.t. probability (p)
        dp = d ./ (p .* log10(p) * log(10));

        % Partial derivative w.r.t. attenuation coefficient (alpha)
        da = -d ./ alpha;
    end
end
