function [x, dp] = logit(p)
    %logit: Computes the logit (log-odds) transformation of a probability p.
    %       Also computes the derivative respect to p if requested.
    %
    % Inputs:
    %     p   - Probability value(s) in the range (0, 1)
    %
    % Outputs:
    %     x   - Logit (log-odds) value(s) for input p
    %     dp  - Derivative(s) of the logit with respect to p
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Compute the logit transformation
    x = log(p ./ (1 - p));

    % Compute the derivative if requested
    if nargout > 1
        dp = 1 ./ (p .* (1 - p));
    end

end
