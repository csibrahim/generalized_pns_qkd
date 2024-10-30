function [p_i, d_pe] = p_bs(i, x, a, b, pe)
    %p_bs: This function computes the probability of a photon exiting path i (0 or 1)
    %      in a beam splitter, based on the bit choice x, basis choices a and b, 
    %      and the beam splitter misalignment probability pe.
    %
    % Inputs:
    %   i  - Path index (0 or 1)
    %   x  - Binary input, Alice's bit choice for the pulse (0 or 1)
    %   a  - Binary input, Alice's basis choice (0 or 1)
    %   b  - Binary input, Bob's basis choice (0 or 1)
    %   pe - Beam splitter misalignment parameter (between 0 and 1)
    %
    % Outputs:
    %   p_i   - Probability that the photon exits path i
    %   d_pe  - (Optional) Derivative of the probability w.r.t. pe
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Compute the angle corresponding to the misalignment parameter
    pe_hat = asin(sqrt(pe));
    
    % Compute the angle theta based on pulse parameters and pe_hat
    theta = (pi/2)*(x + i) - (pi/4)*(a - b) + pe_hat;
    
    % Compute the probability of exiting path i
    p_i = cos(theta).^2;
    
    % If requested, compute the derivative w.r.t. pe
    if nargout > 1
        d_pe = -sin(2 * theta) ./ (2 * sqrt(pe .* (1 - pe)));
    end
end