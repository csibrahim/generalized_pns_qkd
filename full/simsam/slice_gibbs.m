function [samples, logPDFs] = slice_gibbs(logpdf, x0, N, print, varargin)
    %slice_gibbs: Performs Gibbs sampling using a slice sampling method for each dimension.
    %             The method iterates over each dimension of the sample, performing slice sampling 
    %             to explore the distribution defined by the log-pdf function.
    %
    % Inputs:
    %     logpdf   - Function handle that computes the log-pdf at a given point
    %     x0       - Initial sample point (row vector)
    %     N        - Number of samples to draw
    %     print    - Boolean flag for printing progress (optional, default: false)
    %     varargin - Optional parameters:
    %                   'width'   - Step size for slice sampling (scalar or vector, default: 2)
    %                   'maxIter' - Maximum number of steps for stepping out (default: 100)
    %
    % Outputs:
    %     samples  - NxD array of samples drawn from the distribution (N samples, D dimensions)
    %     logPDFs  - Nx1 array of log-pdf values evaluated at the sample points
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    if nargin < 4
        print = false;
    end

    % Parse optional input arguments
    p = inputParser;
    addParameter(p, 'width', 2);
    addParameter(p, 'maxIter', 100);
    parse(p, varargin{:});

    width = p.Results.width;
    maxIter = p.Results.maxIter;

    % Initialize starting point and determine dimensions
    dim = length(x0); % Number of dimensions

    % Handle 'width' parameter: either a scalar or vector with dimension matching x0
    if isscalar(width)
        w = repmat(width, 1, dim);
    elseif isvector(width) && length(width) == dim
        w = width(:)'; % Ensure row vector
    else
        error('Parameter ''width'' must be either a scalar or a vector with length equal to the number of dimensions.');
    end

    % Preallocate arrays for samples and log-pdf values
    samples = zeros(N, dim);
    logPDFs = zeros(N, 1);

    % Evaluate logpdf at the initial point
    logpdf_current = logpdf(x0);

    % Store the initial sample and log-pdf value
    samples(1, :) = x0;
    logPDFs(1) = logpdf_current;

    % Main Sampling Loop
    x_current = x0;
    for i = 2:N
        % Iterate over each dimension
        for d = 1:dim
            % Current value of the d-th dimension
            x_d_current = x_current(d);

            % Sample vertical level: log(y) = log(f(x)) + log(u), u ~ Uniform(0,1)
            log_u = logpdf_current + log(rand());

            % Stepping out to find the slice interval [L, R] for dimension d
            u = rand(); % Uniform random for interval start
            L = x_d_current - w(d) * u; % Left bound
            R = L + w(d);               % Right bound

            % Step out to the left
            j = 0;
            while j < maxIter && logpdf([x_current(1:d-1), L, x_current(d+1:end)]) > log_u
                L = L - w(d); % Expand left
                j = j + 1;
            end

            % Step out to the right
            j = 0;
            while j < maxIter && logpdf([x_current(1:d-1), R, x_current(d+1:end)]) > log_u
                R = R + w(d); % Expand right
                j = j + 1;
            end

            % Shrinkage loop: find a new point within [L, R]
            while true
                % Sample uniformly from [L, R]
                x_new_d = L + (R - L) * rand();

                % Update the proposed sample for the d-th dimension
                x_proposed = x_current;
                x_proposed(d) = x_new_d;

                % Evaluate logpdf at the proposed point
                log_f_new = logpdf(x_proposed);

                if log_f_new >= log_u
                    % Accept the new sample
                    x_current(d) = x_new_d;
                    logpdf_current = log_f_new;
                    break;
                elseif x_new_d < x_current(d)
                    % Shrink the interval from the left
                    L = x_new_d;
                else
                    % Shrink the interval from the right
                    R = x_new_d;
                end
            end
        end

        % Store the updated sample and logpdf
        samples(i, :) = x_current;
        logPDFs(i) = logpdf_current;

        % Display progress if requested
        if print
            progress(i - 1, N - 1, 'sampling');
        end
    end
end