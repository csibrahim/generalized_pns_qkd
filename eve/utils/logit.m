function [x, dp] = logit(p)
    %logit: This function computes the logit (log-odds) transformation of a 
    %        probability p. It can also compute the derivative of the logit 
    %        transformation with respect to p if requested.
    %
    % Inputs:
    %     p   - Probability value(s) in the range (0, 1)
    %
    % Outputs:
    %     x   - Logit (log-odds) values corresponding to the input p
    %     dp  - Derivative of the logit with respect to p
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Compute the logit transformation
    x = log(p ./ (1 - p));

    % If derivatives are requested, compute the derivative
    if nargout > 1
        % Derivative of logit with respect to p
        dp = 1 ./ (p .* (1 - p));
    end

end