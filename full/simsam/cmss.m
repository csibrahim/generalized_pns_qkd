function [samples, logPDFs] = cmss(logpdf, x0, N, varargin)
    %cmss: Covariance-Matching Slice Sampler (CMSS) implementation based on 
    %       the work of Madeleine Thompson and Radford Neal, "Covariance-Adaptive Slice Sampling".
    %       This MATLAB code is a translation of the original R code provided by the authors,
    %       available at https://glizen.com/radfordneal/cass.software.html.
    %
    % Inputs:
    %     logpdf   - Function handle for the log of the target probability density function (PDF)
    %     x0       - Initial sample point
    %     N        - Number of samples to draw
    %     varargin - Optional name-value pair arguments:
    %                   'print'   - Print progress (default: false)
    %                   'sigma_c' - Initial scaling factor for the covariance matrix (default: 2)
    %                   'theta'   - Update rate for the covariance matrix (default: 1)
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
    addParameter(p, 'sigma_c', 2);
    addParameter(p, 'theta', 1);

    % Parse input arguments
    parse(p, varargin{:});

    % Assign parsed values to variables
    print = p.Results.print;
    sigma_c = p.Results.sigma_c;
    theta = p.Results.theta;

    epsilon = 1e-10; % Small constant to avoid division by zero
    
    d = length(x0);  % Dimensionality of the problem
    
    samples = zeros(N, d);  % Preallocate sample matrix
    logPDFs = zeros(N,1);   % Preallocate log-pdf vector

    samples(1, :) = x0;     % Set the first sample
    logPDFs(1) = logpdf(x0); % Evaluate log-pdf at the initial point

    % Main Sampling Loop
    for i = 2:N
        x0 = samples(i-1, :)';  % Previous sample
        current_logpdf = logPDFs(i-1);  % Current log-pdf value
        
        e = exprnd(1);  % Exponentially distributed random variable
        logpdf_y = current_logpdf - e;  % Slice threshold
        
        R = (eye(d) + epsilon * eye(d)) / sigma_c;  % Initial covariance matrix
        F = (eye(d) + epsilon * eye(d)) / sigma_c;  % Initial scaling factor

        cbar_star = zeros(d, 1);  % Un-normalized crumb mean

        % Slice Sampling Loop
        while true
            z = randn(d, 1);  % Standard normal random vector
            c = x0 + F \ z;  % Proposed crumb
            cbar_star = cbar_star + F' * F * c;  % Update the un-normalized crumb mean
            cbar = R \ (R' \ cbar_star);  % Compute the normalized crumb mean

            % Propose new sample
            z = randn(d, 1);
            x = cbar + R \ z;

            % Evaluate log-pdf at the new sample
            [logpdf_x, G] = logpdf(x');

            if logpdf_x >= logpdf_y  % If accepted
                break;
            end

            % Normalize gradient
            G = G';
            norm_G = norm(G);

            if norm_G < epsilon  % Break if gradient is too small
                break;
            end

            g = G / (norm_G + epsilon);

            delta = norm(x - c);  % Step size

            u = x + delta * g;  % Update sample point

            % Evaluate log-pdf at the updated point
            logpdf_u = logpdf(u');

            kappa = -2 * delta^-2 * (logpdf_u - logpdf_x - delta * norm_G) + epsilon;

            % Update the slice region
            lxu = 0.5 * (norm_G^2 / (kappa + epsilon)) + logpdf_x;

            current_logpdf = max(current_logpdf, lxu);

            % Update covariance matrix
            sigma2 = (2 / 3) * ((current_logpdf - logpdf_y) / (kappa + epsilon));
            Rg = R * g; 
            alpha = max(0, (1 / (sigma2 + epsilon)) - (1 + theta) * (Rg' * Rg + epsilon));
            F = cholupdate(sqrt(theta) * R, sqrt(alpha) * g, '+');
            R = cholupdate(sqrt(theta + 1) * R, sqrt(alpha) * g, '+');
        end

        % Store the accepted sample and log-pdf value
        samples(i, :) = x';
        logPDFs(i) = logpdf_x;

        if print
            progress(i-1, N-1, 'sampling');
        end
    end
end