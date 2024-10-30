function [samples, logPDFs] = srss(logpdf, x0, N, varargin)
    %srss: Shrinking-Rank Slice Sampler (SRSS) implementation based on
    %       the work of Thompson and Neal (2010), "Covariance-Adaptive Slice Sampling".
    %       This MATLAB code is an adaptation of the original R code by the authors,
    %       available at https://glizen.com/radfordneal/cass.software.html.
    %
    % Inputs:
    %     logpdf    - Function handle for the log of the target probability density function (PDF)
    %     x0        - Initial sample point
    %     N         - Number of samples to draw
    %     varargin  - Optional name-value pair arguments:
    %                       'print'         - Print progress (default: false)
    %                       'sigma_c'       - Initial crumb scale (default: 1)
    %                       'downscale'     - Shrinkage factor for the crumb scale (default: 0.9)
    %                       'min_dimension' - Minimum dimensionality for projection (default: 1)
    %
    % Outputs:
    %     samples  - Matrix of sampled points (N x d)
    %     logPDFs  - Vector of log-pdf values at each sampled point (N x 1)
    %
    % Reference:
    %     Thompson and Neal (2010), "Covariance-Adaptive Slice Sampling"
    %     (https://arxiv.org/abs/1003.3201).
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Set up input parser with default values
    p = inputParser;
    addParameter(p, 'print', false);
    addParameter(p, 'sigma_c', 1);
    addParameter(p, 'downscale', 0.9);
    addParameter(p, 'min_dimension', 1);

    % Parse input arguments
    parse(p, varargin{:});

    % Assign parsed values to variables
    print = p.Results.print;
    sigma_c = p.Results.sigma_c;
    downscale = p.Results.downscale;
    min_dimension = p.Results.min_dimension;

    % Projection function
    P = @(M, a) a - (a * M') * M;

    d = length(x0);  % Dimensionality of the problem
    samples = NaN(N, d);  % Preallocate sample matrix
    logPDFs = NaN(N, 1);  % Preallocate log-pdf vector

    samples(1, :) = x0(:)';  % Set the first sample
    logPDFs(1) = logpdf(x0);  % Evaluate log-pdf at the initial point

    % Main Sampling Loop
    for i = 2:N
        x_prev = samples(i - 1, :);  % Previous sample
        logpdf_x0 = logpdf(x_prev);  % Log-pdf at previous sample
        logpdf_y = logpdf_x0 - exprnd(1);  % Slice threshold

        J = zeros(0, d);  % Initialize projection matrix
        k = 0;
        cbar_star = zeros(1, d);  % Un-normalized mean of crumbs
        Lambda = 0;  % Proposal precision

        % Slice Sampling Loop
        while true
            k = k + 1;
            sigma_ck = sigma_c * downscale^(k - 1);  % Crumb standard deviation
            crumb = sigma_ck * randn(1, d);  % Generate crumb
            Lambda = Lambda + sigma_ck^(-2);  % Update precision
            cbar_star = cbar_star + crumb / sigma_ck^2;  % Update crumb mean

            % Propose new sample
            proposal = cbar_star / Lambda + randn(1, d) / sqrt(Lambda);
            x = x_prev + P(J, proposal);

            % Evaluate log-pdf at the proposed point
            logpdf_x = logpdf(x);

            if logpdf_x >= logpdf_y  % Accept the proposal
                break;
            end

            % Update projection matrix if possible
            if (size(J, 1) < d - min_dimension)
                [~, d_logpdf_x] = logpdf(x);  % Gradient at proposed point
                gx = d_logpdf_x;  % Gradient vector
                gstar = P(J, gx);  
                norm_gstar = norm(gstar);
                norm_gx = norm(gx);

                if (norm_gstar == 0 || norm_gx == 0)
                    cosine_similarity = 0;
                else
                    cosine_similarity = (gstar * gx') / (norm_gstar * norm_gx);
                end

                if (cosine_similarity > 0.5)
                    new_column = P(J, gstar / norm_gstar);  % Normalize and add to J
                    J = [J; new_column];  % Append as new row
                end
            end
        end

        samples(i, :) = x;
        logPDFs(i) = logpdf_x;

        % Update progress if requested
        if (print)
            progress(i-1, N-1, 'sampling');
        end
    end
end