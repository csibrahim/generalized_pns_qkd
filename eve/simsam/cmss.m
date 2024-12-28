function [samples, logPDFs] = cmss(logpdf, x0, N, varargin)
    %cmss: Covariance-Matching Slice Sampler (CMSS) implementation based on
    %      the work of Madeleine Thompson and Radford Neal, "Covariance-Adaptive 
    %      Slice Sampling". This MATLAB implementation translates the original
    %      R code provided by the authors, available at:
    %      https://glizen.com/radfordneal/cass.software.html.
    %
    % Inputs:
    %     logpdf   - Function handle for the log of the target probability density
    %     x0       - Initial sample point (vector of size d)
    %     N        - Number of samples to generate
    %     varargin - Optional name-value pair arguments:
    %                   'print'   - Display sampling progress (default: false)
    %                   'sigma_c' - Initial scaling factor for the covariance matrix (default: 2)
    %                   'theta'   - Update rate for adapting the covariance matrix (default: 1)
    %
    % Outputs:
    %     samples  - N x d matrix of sampled points
    %     logPDFs  - N x 1 vector of log-pdf values at each sampled point
    %
    % Reference:
    %     Thompson and Neal (2010), "Covariance-Adaptive Slice Sampling"
    %     (https://arxiv.org/abs/1003.3201).
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Set up input parser with default values
    p = inputParser;
    addParameter(p, 'print', false);    % Progress display flag
    addParameter(p, 'sigma_c', 2);     % Initial covariance scaling
    addParameter(p, 'theta', 1);       % Covariance update rate

    % Parse input arguments
    parse(p, varargin{:});

    % Assign parsed values to variables
    print = p.Results.print;          % Display sampling progress
    sigma_c = p.Results.sigma_c;      % Initial covariance scaling factor
    theta = p.Results.theta;          % Covariance matrix update rate

    epsilon = 1e-10; % Small constant to prevent division by zero
    
    d = length(x0);  % Dimensionality of the problem
    
    % Preallocate storage for samples and log-pdf values
    samples = zeros(N, d);  % Sample storage matrix
    logPDFs = zeros(N, 1);  % Log-pdf storage vector

    % Initialize the first sample and log-pdf value
    samples(1, :) = x0;      
    logPDFs(1) = logpdf(x0);

    % Main Sampling Loop
    for i = 2:N
        x0 = samples(i-1, :)';  % Previous sample
        current_logpdf = logPDFs(i-1);  % Log-pdf value at previous sample
        
        e = exprnd(1);  % Exponentially distributed random variable for slice threshold
        logpdf_y = current_logpdf - e;  % Slice threshold for the current iteration
        
        % Initialize covariance matrix and scaling factor
        R = (eye(d) + epsilon * eye(d)) / sigma_c;  
        F = (eye(d) + epsilon * eye(d)) / sigma_c;  

        cbar_star = zeros(d, 1);  % Accumulator for un-normalized crumb mean

        % Slice Sampling Loop
        while true
            z = randn(d, 1);  % Draw standard normal random vector
            c = x0 + F \ z;   % Proposed crumb
            
            % Update the un-normalized crumb mean
            cbar_star = cbar_star + F' * F * c;  
            cbar = R \ (R' \ cbar_star);  % Compute the normalized crumb mean

            % Propose a new sample
            z = randn(d, 1);  % Generate a new normal random vector
            x = cbar + R \ z; % Proposed new sample

            % Evaluate log-pdf at the proposed sample
            [fx, gx] = logpdf(x');

            if fx >= logpdf_y  % Accept sample if within the slice
                break;
            end

            % Normalize gradient of log-pdf
            norm_gx = norm(gx');

            if norm_gx < epsilon  % Stop if gradient norm is too small
                break;
            end

            g_norm = gx' / (norm_gx + epsilon);  % Normalized gradient

            % Calculate step size and update sample point
            delta = norm(x - c);  
            u = x + delta * g_norm;    

            % Evaluate log-pdf at the updated point
            logpdf_u = logpdf(u');

            % Compute the curvature parameter
            kappa = -2 * delta^-2 * (logpdf_u - fx - delta * norm_gx) + epsilon;

            % Update slice region
            lxu = 0.5 * (norm_gx^2 / (kappa + epsilon)) + fx;
            current_logpdf = max(current_logpdf, lxu);

            % Update covariance matrix with gradient-based adaptation
            sigma2 = (2 / 3) * ((current_logpdf - logpdf_y) / (kappa + epsilon));
            Rg = R * g_norm; 
            alpha = max(0, (1 / (sigma2 + epsilon)) - (1 + theta) * (Rg' * Rg + epsilon));
            F = cholupdate(sqrt(theta) * R, sqrt(alpha) * g_norm, '+');      % Update scaling
            R = cholupdate(sqrt(theta + 1) * R, sqrt(alpha) * g_norm, '+');  % Update covariance
        end

        % Store the accepted sample and log-pdf value
        samples(i, :) = x';
        logPDFs(i) = fx;

        % Display progress if enabled
        if print
            progress(i-1, N-1, 'sampling');
        end
    end
end
