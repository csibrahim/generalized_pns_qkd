function [samples, logPDFs] = slice_oneshot(logpdf, x0, N, varargin)
    %slice_oneshot: Generates samples from a target distribution using slice sampling.
    %
    % This function is based on MATLAB's slice sampling algorithm, with additional 
    % functionality to retrieve the log-probabilities (`logPDFs`) of sampled points 
    % and optionally display progress throughout the sampling process. The log-probabilities 
    % are calculated using the user-provided `logpdf` function handle.
    %
    % Inputs:
    %   logpdf   - Function handle for the log of the target probability density function (PDF)
    %   x0       - Initial sample point (row vector, 1 x d)
    %   N        - Number of samples to draw
    %   varargin - Optional name-value pair arguments:
    %                'print'   - Boolean to display progress (default: false)
    %                'width'   - Step size for slice expansion (default: 2)
    %                'maxIter' - Maximum number of iterations for stepping out/shrinking (default: none)
    %
    % Outputs:
    %   samples  - N x d matrix of sampled points, where N is the number of samples and d is the dimension
    %   logPDFs  - N x 1 vector of log-probabilities of the sampled points, evaluated at each sample
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Parse optional input arguments
    p = inputParser;
    addParameter(p, 'print', false);
    addParameter(p, 'width', 2);
    addParameter(p, 'maxIter', Inf);  % Default for maxIter is Inf if not specified
    parse(p, varargin{:});

    print = p.Results.print;
    width = p.Results.width;
    maxIter = p.Results.maxIter;

    dim = size(x0, 2);           % Dimension of the distribution
    samples = zeros(N, dim);     % Placeholder for the random sample sequence
    logPDFs = zeros(N, 1);       % Placeholder for the logpdf values
    
    e = exprnd(1, N, 1);         % Vertical position of the slice
    RW = rand(N, dim);           % Randomized width factors
    RD = rand(N, dim);           % Random points within the slice interval
    
    % Function to check if a point is within the slice (i.e., above threshold)
    inside = @(x, th) (logpdf(x) > th); 
    
    % Sampling loop
    for i = 1:N
        % Define the slice threshold using a vertical level drawn from (0, f(x0))
        z = logpdf(x0) - e(i);

        % Random initial interval around x0 for stepping-out procedure
        r = width .* RW(i, :); % randomized width/step size
        xl = x0 - r; 
        xr = xl + width; 
        iter = 0;

        % Perform stepping-out procedure if univariate
        if dim == 1
            % Step out to the left
            while inside(xl, z) && iter < maxIter
                xl = xl - width;
                iter = iter + 1;
            end
            if iter >= maxIter || any(xl < -sqrt(realmax))
                error("Stepping out took too many iterations.");
            end

            % Step out to the right
            iter = 0;
            while inside(xr, z) && iter < maxIter
                xr = xr + width;
                iter = iter + 1;
            end
            if iter >= maxIter || any(xr > sqrt(realmax))
                error("Stepping out took too many iterations.");
            end
        end

        % Draw a point uniformly from the interval [xl, xr]
        xp = RD(i, :) .* (xr - xl) + xl;

        % Shrink the interval if the point is outside the slice
        iter = 0;
        while ~inside(xp, z) && iter < maxIter
            rshrink = (xp > x0);
            xr(rshrink) = xp(rshrink);
            lshrink = ~rshrink;
            xl(lshrink) = xp(lshrink);
            xp = rand(1, dim) .* (xr - xl) + xl;
            iter = iter + 1;
        end

        if iter >= maxIter
            error("Shrinkage procedure took too many iterations.");
        end

        % Update current sample and store the log-pdf value
        x0 = xp;
        samples(i, :) = x0;
        logPDFs(i) = logpdf(x0);

        % Display progress if 'print' option is enabled
        if print
            progress(i, N, 'Sampling progress');
        end
    end
end