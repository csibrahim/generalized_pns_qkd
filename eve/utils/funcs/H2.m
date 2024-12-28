function entropy = H2(p)
    %H2: Binary Shannon entropy for probability p.
    %
    % Inputs:
    %     p        - Probability value(s) (scalar, vector, or matrix)
    %
    % Outputs:
    %     entropy  - Binary Shannon entropy value(s) for input p
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Compute binary Shannon entropy, handling p = 0 or p = 1 edge cases
    entropy = -p .* log2(p + (p == 0)) - (1 - p) .* log2(1 - p + (p == 1));

end