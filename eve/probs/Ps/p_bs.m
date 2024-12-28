function [p_i, d_pe] = p_bs(i, x, a, b, pe)
    %p_bs: Computes the probability of a photon exiting path (i), 0 or 1, 
    %      in a beam splitter based on Alice's bit (x), basis choices 
    %      (a, b), and the beam splitter misalignment probability (pe). 
    %      Optionally computes the derivative of the probability w.r.t. pe.
    %
    % Inputs:
    %   i  - Binary input, Path index (0 or 1)
    %   x  - Binary input, Alice's bit choice (0 or 1)
    %   a  - Binary input, Alice's basis choice (0 or 1)
    %   b  - Binary input, Bob's basis choice (0 or 1)
    %   pe - Beam splitter misalignment probability (0 ≤ pe ≤ 1)
    %
    % Outputs:
    %   p_i   - Probability that the photon exits path i
    %   d_pe  - Derivative of p_i w.r.t. pe
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Compute the angle corresponding to the misalignment probability pe
    pe_hat = asin(sqrt(pe));

    % Compute the angle theta based on the input parameters and pe_hat
    theta = (pi / 2) * (x + i) - (pi / 4) * (a - b) + pe_hat;

    % Compute the probability of exiting path i
    p_i = cos(theta).^2;

    % Compute the derivative of the p_i w.r.t. pe if requested
    if nargout > 1
        d_pe = -sin(2 * theta) ./ (2 * sqrt(pe .* (1 - pe)));
    end
end
