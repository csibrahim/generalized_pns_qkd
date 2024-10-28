function entropy = H2(p)
    %H2: This function computes the binary Shannon entropy H(p) for a given 
    %     probability p. It handles edge cases where p = 0 or p = 1 to avoid 
    %     undefined logarithms.
    %
    % Inputs:
    %     p        - Probability value(s) (can be scalar, vector, or matrix)
    %
    % Outputs:
    %     entropy  - Binary Shannon entropy value(s) corresponding to input p
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Compute the binary Shannon entropy, handling edge cases where p is 0 or 1
    entropy = -p .* log2(p + (p == 0)) - (1 - p) .* log2(1 - p + (p == 1));

end