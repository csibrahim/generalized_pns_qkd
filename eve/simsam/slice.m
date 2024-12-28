function [samples, logPDFs] = slice(logpdf, x0, N, varargin)
    %slice: Generates samples using a modified slice sampling algorithm.
    %       This implementation is based on MATLAB's built-in `slicesample`
    %       function, with added flexibility for controlling the output. It
    %       calculates log-probabilities (`logPDFs`) of sampled points and
    %       allows for progress display during sampling.
    %
    % Inputs:
    %   logpdf   - Handle for the log of the target probability density function
    %   x0       - Initial sample point (row vector, 1 x d)
    %   N        - Number of samples to generate
    %   varargin - Optional name-value pairs:
    %                'print'   - Display progress during sampling (default: false)
    %                'width'   - Step size for slice expansion (default: 2)
    %                'maxIter' - Max iterations for stepping out/shrinking (default: 200)
    %
    % Outputs:
    %   samples  - N x d matrix of sampled points
    %   logPDFs  - N x 1 vector of log-probabilities for sampled points
    %
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Parse optional input arguments
    p = inputParser;
    addParameter(p, 'print', false);
    addParameter(p, 'width', 2);
    addParameter(p, 'maxIter', 200);
    parse(p, varargin{:});

    % Extract parsed parameters
    print = p.Results.print;
    width = p.Results.width;
    maxIter = p.Results.maxIter;

    % Initializations
    d = size(x0, 2);           % Dimensionality of the distribution
    samples = zeros(N, d);     % Preallocate samples matrix
    logPDFs = zeros(N, 1);     % Preallocate log-PDFs vector
    
    e = exprnd(1, N, 1);       % Slice heights
    RW = rand(N, d);           % Randomized width multipliers
    RD = rand(N, d);           % Random points for slice interval
    
    % Function to check if a point is within the slice
    inside = @(x, th) (logpdf(x) > th); 
    
    % Sampling loop
    for i = 1:N
        % Determine slice threshold
        z = logpdf(x0) - e(i);

        % Initialize slice interval around x0
        r = width .* RW(i, :); % Randomized step size
        xl = x0 - r; 
        xr = xl + width; 
        iter = 0;

        % Perform stepping-out procedure for univariate case
        if d == 1
            % Step out to the left
            while inside(xl, z) && iter < maxIter
                xl = xl - width;
                iter = iter + 1;
            end
            if iter >= maxIter
                error("Stepping out exceeded maximum iterations.");
            end

            % Step out to the right
            iter = 0;
            while inside(xr, z) && iter < maxIter
                xr = xr + width;
                iter = iter + 1;
            end
            if iter >= maxIter
                error("Stepping out exceeded maximum iterations.");
            end
        end

        % Draw a point uniformly from [xl, xr]
        xp = RD(i, :) .* (xr - xl) + xl;

        % Shrink interval if point is outside the slice
        iter = 0;
        while ~inside(xp, z) && iter < maxIter
            rshrink = (xp > x0);  % Determine shrink direction
            xr(rshrink) = xp(rshrink);
            lshrink = ~rshrink;
            xl(lshrink) = xp(lshrink);
            xp = rand(1, d) .* (xr - xl) + xl;
            iter = iter + 1;
        end

        if iter >= maxIter
            error("Shrinkage exceeded maximum iterations.");
        end

        % Update current sample and store log-pdf value
        x0 = xp;
        samples(i, :) = x0;
        logPDFs(i) = logpdf(x0);

        % Display progress if requested
        if print
            progress(i, N, 'Sampling progress');
        end
    end
end
