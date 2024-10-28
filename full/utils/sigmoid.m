function [p, dx] = sigmoid(x)
    %sigmoid: This function computes the sigmoid (logistic) function of the input x,
    %          and optionally returns its derivative with respect to x.
    %
    % Inputs:
    %     x   - Input value(s) (can be scalar, vector, or matrix)
    %
    % Outputs:
    %     p   - Sigmoid function values corresponding to the input x
    %     dx  - Derivative of the sigmoid function with respect to x
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Compute the sigmoid function
    p = 1 ./ (1 + exp(-x));

    % Compute the derivative if requested
    if nargout > 1
        % Derivative of the sigmoid function with respect to x
        dx = p .* (1 - p);
    end

end