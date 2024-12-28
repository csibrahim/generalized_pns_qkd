function plotPMFs(N, p, CIp, data, FontSize, label)
    %plotPMFs: Creates a bar plot comparing observed data with theoretical 
    %          binomial probabilities, including confidence intervals.
    %
    % Inputs:
    %   N         - Total trials for the binomial distribution
    %   p         - Success probability for each trial
    %   CIp       - Confidence level (e.g., 0.95 for 95% CI)
    %   data      - Observed data values
    %   FontSize  - Font size for labels (default: 16)
    %   label     - X-axis label (default: none)
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).
   
    % Set default values for optional inputs
    if (nargin < 5),  FontSize = 16; end
    if (nargin < 6),  label = []; end

    % Adjust font size for tick labels
    tickFontSize = 0.5 * FontSize;

    % Compute confidence interval bounds
    epsilon = (1 - CIp) / 2;
    lowers = binoinv(epsilon, N, p);
    uppers = binoinv(1 - epsilon, N, p);

    % Calculate marker size and line width based on plot dimensions
    current_units = get(gca, 'Units');
    set(gca, 'Units', 'points');
    position = get(gca, "Position");

    marker_size = position(3) / 20;
    line_width = position(3)/300;

    set(gca, 'Units', current_units);

    % Define range for plot limits
    from = min([data;lowers]);
    to = max([data;uppers]);

    % Compute theoretical PMF
    n = from:to;
    pmf = binopdf(n, N, p);

    % Compute histogram counts
    [counts, edges] = histcounts(data, 'Normalization', 'pdf');

    % Compute bin centers for plotting
    binCenters = edges(1:end-1) + diff(edges) / 2;
    
    % Determain plotting upper-limit
    top = max([counts(:); pmf(:)]) * 1.1;

    hold on;

    % Plot confidence interval region
    area([lowers-0.5, uppers+0.5], [1, 1] * top, ...
         'FaceColor', [1, 1, 1] * 0.8, ...
         'EdgeColor', 'k');

    % Plot the histogram using bar
    h_hist = bar(binCenters, counts, ...
                 'FaceColor', 'w', ...
                 'FaceAlpha', 1, ...
                 'LineWidth', line_width...
                 );

    % Plot theoretical PMF as markers
    plot(n, pmf, 'o', ...
         'MarkerSize', marker_size, ...
         'MarkerFaceColor', 'k', ...
         'MarkerEdgeColor', 'none');
    

    % Set axis limits
    ylim([0 top]);
    xlim([from-0.5 to+0.5]);

    % Format tick labels and adjust font sizes
    set(gca, ...
        'TickLabelInterpreter', 'latex', ...
        'FontSize', tickFontSize...
        );

    ax = gca;
    ax.YAxis.FontSize = tickFontSize;
    ax.XAxis.FontSize = tickFontSize;

    % Add x-axis label
    xlabel(label, ...
           'Interpreter', 'latex', ...
           'FontSize', FontSize ...
           );

    % Adjust x-ticks to show integers only
    currentXTicks = xticks;
    
    if (max(currentXTicks) <= 5)
        % If less than 5, display ticks for all integers
        xticks(0:max(currentXTicks));
    else
        % Retain only the integer ticks
        integerXTicks = currentXTicks(currentXTicks == floor(currentXTicks));
        xticks(integerXTicks);
    end

end