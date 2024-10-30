function [samples, logPDFs] = sample(method, theta0, Ns, Nb, C, thetaA, thetaB, thetaP, varargin)
    %sample: Generates samples using a specified sampling method, with optional
    %         burn-in and display of progress. Allows interruption during sampling.
    %
    % Inputs:
    %     method    - Sampler options:
    %                    'cmss'  - Covariance-Matching Slice Sampler
    %                    'srss'  - Shrinking-Rank Slice Sampler
    %                    'slice' - Slice sampling (within Gibbs)
    %     theta0    - Initial parameter vector
    %     Ns        - Number of sampling iterations
    %     Nb        - Number of burn-in iterations
    %     C         - Counts or observations for likelihood calculation
    %     thetaA    - Cell array of system parameters specific to the channel between 
    %                 Alice and Bob, containing [lambdas, alpha, dAB].
    %     thetaB    - Cell array of system parameters for Bob's detectors, containing 
    %                 [pa0, pa1, pc0, pc1, pd0, pd1, pe].
    %     thetaP    - Prior parameters as a cell array {alphas, betas, ub, lb}
    %     varargin  - Optional name-value pair arguments:
    %                       'chunkSize'    - Number of samples processed per chunk (default: 10)
    %                       'ground_truth' - Ground truth parameter values (default: [])
    %                       'display'      - Boolean for progress display (default: true)
    %                       'labels'       - Labels for plotting (default: [])
    %                       'FontSize'     - Font size for labels and text (default: 30)
    %                       'FigureWidth'  - Width of the figure in points (default: 1035)
    %
    % Outputs:
    %     samples  - Sampled values for each parameter after burn-in
    %     logPDFs  - Log-probability of each sampled state
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).
    
    % Set up input parser with default values
    p = inputParser;
    addParameter(p, 'chunkSize', 10);
    addParameter(p, 'ground_truth', []);
    addParameter(p, 'display', true);
    addParameter(p, 'labels', []);
    addParameter(p, 'FontSize', 30);
    addParameter(p, 'FigureWidth', 1035);

    % Parse input arguments
    parse(p, varargin{:});

    % Assign parsed values to variables
    chunkSize = p.Results.chunkSize;
    ground_truth = p.Results.ground_truth;
    display = p.Results.display;
    labels = p.Results.labels;
    FontSize = p.Results.FontSize;
    FigureWidth = p.Results.FigureWidth;

    % Total iterations, including burn-in
    N = Ns + Nb;

    % Convert initial parameters to their transformed space using prior parameters
    phi0 = Phi(theta0, thetaP);

    % Define logpdf as a function of transformed parameters for sampling
    logpdf = @(phiE) logPDF(C, thetaA, thetaB, phiE, thetaP);
    
    % Initialize storage for samples and log probabilities
    dim = length(phi0);
    samples = zeros(N, dim);
    logPDFs = zeros(N, 1);

    chunks = ceil(N / chunkSize); % Number of chunks to process
    chunk = phi0;                 % Initialize chunk with starting parameter values

    shouldBreak = false; % Flag to handle interruption by user

    % Create figure for displaying progress, if enabled
    if display
        f = prepFigure(FigureWidth, FontSize);
        set(f, 'KeyPressFcn', @keyPressCallback);
    end

    % Start sampling process
    tic;
    for i = 1:chunks
        % Exit loop if interrupted
        if (shouldBreak)
            disp('Loop interrupted by user command. Exiting gracefully.');
            break;
        end

        % Define chunk range
        first = 1 + (i - 1) * chunkSize;
        last = min(i * chunkSize, N);

        % Perform sampling based on selected method
        switch method
            case 'cmss'
                [chunk, logPDFs(first:last)] = cmss(logpdf, chunk(end, :), chunkSize);
            case 'srss'
                [chunk, logPDFs(first:last)] = srss(logpdf, chunk(end, :), chunkSize);
            case 'slice'
                [chunk, logPDFs(first:last)] = slice_gibbs(logpdf, chunk(end, :), chunkSize);
        end

        % Convert sampled values from transformed space back to original space
        samples(first:last, :) = Theta(chunk, thetaP);

        % Display progress and histograms if enabled
        if display
            time = toc;
            if (time > 1 / 10) % Update display every 0.1 seconds
                sampling_percentage = max((i * chunkSize - Nb) / Ns, 0);
                sampling_percentage = round(1e4 * sampling_percentage) / 1e2;

                burnin_percentage = min(i * chunkSize / Nb, 1);
                burnin_percentage = round(1e4 * burnin_percentage) / 1e2;

                % Display progress on figure title
                set(f, 'Name', ['burn-in ', num2str(burnin_percentage), '% - sampling ', num2str(sampling_percentage), '%'], 'NumberTitle', 'off');

                % Display histograms for samples (ignoring burn-in samples gradually)
                start = min(max((i - 1) * chunkSize - Nb + 1, 1), Nb);
                finish = last;
                displayHistograms(samples(start:finish, :), thetaP{3}, thetaP{4}, 'ground_truth', ground_truth, 'labels', labels, 'newFigure', false);

                % Annotation with stop command instruction
                width = 1;
                height = 0.05;
                
                annotation('textbox', [0, 0, width, height], ...
                           'String', 'Press CTRL+C any time to stop the sampling process', ...
                           'Color', [0.2 0.2 0.2], ...
                           'HorizontalAlignment', 'center', ...
                           'VerticalAlignment', 'top', ...
                           'FontSize', 12, ...
                           'EdgeColor', 'none', ...
                           'Interpreter', 'latex');
                
                drawnow;
            end
        else
            % Update progress in command window if display is disabled
            progress(i, chunks, 'sampling');
        end
    end

    % Callback function to handle keyboard interruption
    function keyPressCallback(~, event)
        if (strcmp(event.Key, 'c') && ismember('control', event.Modifier))
            shouldBreak = true;
            % Truncate samples and log probabilities to the last complete chunk
            samples = samples(1:(last - chunkSize), :);
            logPDFs = logPDFs(1:(last - chunkSize));
        end
    end

    % If interrupted after burn-in, discard them, otherwise return waht is sampled so far
    if (length(samples) > Nb)
        samples = samples(Nb + 1:end, :);
        logPDFs = logPDFs(Nb + 1:end);
    end
end