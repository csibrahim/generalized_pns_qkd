function [p, dx] = sigmoid(x)
    %sigmoid: Computes the sigmoid (logistic) function for the input x.
    %         Optionally returns the derivative with respect to x.
    %
    % Inputs:
    %     x   - Input value(s) (scalar, vector, or matrix)
    %
    % Outputs:
    %     p   - Sigmoid function value(s) for the input x
    %     dx  - Derivative(s) of the sigmoid function with respect to x
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Compute the sigmoid function
    p = 1 ./ (1 + exp(-x));

    % Compute the derivative if requested
    if nargout > 1
        dx = p .* (1 - p);
    end

end
