function [samples, logPDFs] = srss(logpdf, x0, N, varargin)
    %srss: Shrinking-Rank Slice Sampler (SRSS) implementation based on
    %      the method by Thompson and Neal from "Covariance-Adaptive
    %      Slice Sampling". This MATLAB code translates their original
    %      R implementation, available at:
    %      https://glizen.com/radfordneal/cass.software.html.
    %
    % Inputs:
    %     logpdf    - Handle for the log of the target probability density
    %     x0        - Initial sample (vector of size d)
    %     N         - Number of samples to generate
    %     varargin  - Optional name-value pairs:
    %                   'print'         - Display progress (default: false)
    %                   'sigma_c'       - Initial crumb scale (default: 2)
    %                   'downscale'     - Crumb scale reduction factor (default: 0.9)
    %                   'min_dimension' - Minimum projection dimension (default: 1)
    %                   'maxIter'       - Maximum inner loop iterations (default: 200)
    %
    % Outputs:
    %     samples  - Matrix of generated samples (N x d)
    %     logPDFs  - Log-PDF values at each sampled point (N x 1)
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
    addParameter(p, 'maxIter', 200);
    
    % Parse input arguments
    parse(p, varargin{:});
    
    % Assign parsed values to variables
    print = p.Results.print;
    sigma_c = p.Results.sigma_c;
    downscale = p.Results.downscale;
    min_dimension = p.Results.min_dimension;
    maxIter = p.Results.maxIter;
    
    % Dimensionality of the problem
    d = length(x0);
    samples = NaN(N, d);
    logPDFs = NaN(N, 1);
    
    samples(1, :) = x0(:)';
    logPDFs(1) = logpdf(x0);
    
    % Preallocate J and define maximum rows
    max_J_rows = d - min_dimension;
    J = zeros(max_J_rows, d);  % Preallocated J
    
    % Main Sampling Loop
    for i = 2:N
        
        x_prev = samples(i - 1, :); % Retrieve the previous sample
        logpdf_x0 = logPDFs(i - 1); % Retrieve the log-pdf value at the previous sample

        % Draw a random slice threshold using exponential distribution
        logpdf_y = logpdf_x0 - exprnd(1); 
        
        % Initialize values
        k = 0;                   % Crumb iteration counter
        cbar_star = zeros(1, d); % Cumulative crumb-weighted sum (cbar_star)
        Lambda = 0;              % Proposal precision
        accept = false;          % Acceptance flag
        J_index = 0;             % Index for projection matrix J
        
        % Start the inner loop for crumb sampling and proposal generation
        while ~accept && k < maxIter

            % Increment iteration count
            k = k + 1;
            
            % Compute crumb scale for this iteration
            sigma_ck = sigma_c * downscale^(k - 1);

            % Compute the inverse square of the crumb scale
            sigma_ck_inv2 = sigma_ck^(-2);

            % Generate a random crumb from a normal distribution
            crumb = sigma_ck * randn(1, d);

            % Update precision sum
            Lambda = Lambda + sigma_ck_inv2;

            % Update the weighted sum of crumbs
            cbar_star = cbar_star + crumb * sigma_ck_inv2;
            
            % Generate a proposal 
            mu = cbar_star / Lambda;  % Mean of the proposal distribution
            sigma = 1 / sqrt(Lambda); % Standard deviation of the proposal distribution
            
            proposal = mu + sigma * randn(1, d); % Sample from the proposed distribtion
            
            % Project the proposal if there are active directions in J
            if J_index > 0
                % Retrieve the current projection subspace
                J_current = J(1:J_index, :);

                % Project the proposal onto the subspace
                proj = (proposal * J_current') * J_current;

                % Subtract the projection
                proposal_proj = proposal - proj;
            else
                % No projection needed if J is empty
                proposal_proj = proposal;
            end

            % Translate the proposal to the current state
            x = x_prev + proposal_proj;
            
            % Evaluate logpdf at x
            if J_index < max_J_rows
                % Compute both log-pdf and its gradient
                [fx, gx] = logpdf(x);
            else
                % Compute only log-pdf if J is fully populated
                fx = logpdf(x);
            end
            
            % Accept the proposal if it is inside the slice and break
            if fx >= logpdf_y
                accept = true;
                break;
            end
            
            % Update the projection matrix J
            if J_index < max_J_rows
                
                if J_index > 0
                    % Project the gradient onto the current subspace
                    J_current = J(1:J_index, :);
                    proj_gx = (gx * J_current') * J_current;
                    gstar = gx - proj_gx;
                else
                    % Use the gradient if J is empty
                    gstar = gx;
                end
                
                % Compute the norm of the projected gradient
                norm_gstar = norm(gstar);

                % Compute the norm of the gradient
                norm_gx = norm(gx);
                
                if (norm_gstar == 0 || norm_gx == 0)
                    % Handle edge cases for zero-norm gradients
                    cosine_similarity = 0;
                else
                    % Compute the cosine similarity between gradients
                    cosine_similarity = (gstar * gx') / (norm_gstar * norm_gx);
                end
                
                % Update J if the angle less than 60 degrees
                if (cosine_similarity > 0.5)

                    % Normalize the projected gradient
                    gstar_normed = gstar / norm_gstar;
                    
                    if J_index > 0
                        % Project the normalized gradient onto J
                        J_current = J(1:J_index, :);
                        proj_gstar_normed = (gstar_normed * J_current') * J_current;
                        new_row = gstar_normed - proj_gstar_normed;
                    else
                        % No projection needed if J is empty
                        new_row = gstar_normed;
                    end
                    
                    % Add the new direction to the projection matrix
                    J_index = J_index + 1;
                    J(J_index, :) = new_row;
                end
            end
        end
        
        % If no valid proposal is found within maxIter, revert to the previous sample
        if ~accept
            x = x_prev;
            fx = logpdf_x0;
        end
        
        % Store the accepted sample and its logPDF
        samples(i, :) = x;
        logPDFs(i) = fx;
        
        % Display progress if requested
        if print
            progress(i - 1, N - 1, 'sampling');
        end
    end
end