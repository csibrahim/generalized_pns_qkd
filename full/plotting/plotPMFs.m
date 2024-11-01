function plotPMFs(N, p, epsilon, data, rows, row, cols, col, FontSize)
    %plotPMFs: Generates bar plots with error bars by computing theoretical values and confidence intervals
    %          of a binomial distribution internally from N, p, and epsilon.
    %
    % Inputs:
    %   N           - Number of trials (scalar)
    %   p           - Vector of probabilities for the binomial distribution (theoretical probabilities)
    %   epsilon     - Significance level for confidence intervals
    %   data        - Matrix of data values (rows represent samples, columns represent different variables)
    %   rows        - Number of subplot rows (default: 1)
    %   row         - Current row index (default: 1)
    %   cols        - Number of subplot columns (default: number of data columns)
    %   col         - Current column index (default: 1)
    %   FontSize    - Font size for axis labels and ticks (default: 20)
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Set default values for optional inputs
    if (nargin < 5),  rows = 1; end
    if (nargin < 6),  row = 1; end
    if (nargin < 7),  cols = size(data, 2); end
    if (nargin < 8),  col = 1; end
    if (nargin < 9),  FontSize = 20; end

    % Define derived font sizes for ticks and legends
    legendFontSize = FontSize * 0.75;
    tickFontSize = FontSize * 0.5;

    % Get the number of samples (Ns) and number of variables (Nl)
    Nl = size(data,2);

    % Generate labels for lambda parameters
    labels = arrayfun(@(x) sprintf('$\\lambda_{%d}$', x), 1:Nl, 'UniformOutput', false);

    % Compute confidence intervals
    lowers = binoinv(epsilon, N, p);
    uppers = binoinv(1 - epsilon, N, p);
    
    % Loop over each variable to create subplots
    for i = 1:Nl
        subplot(rows, cols, (row - 1) * cols + (col - 1) + i);

        % Capture the current axes position for setting cap size
        current_units = get(gca, 'Units');
        set(gca, 'Units', 'points');
        pos = get(gca, "Position");

        marker_size = pos(3) / 20;
        line_width = pos(3)/300;

        set(gca, 'Units', current_units);

        % Define error bar sizes relative to plot dimensions
        from = min([data(:,i);lowers(i)]);
        to = max([data(:,i);uppers(i)]);

        area([lowers(i)-0.5, uppers(i)+0.5], [1.1, 1.1], 'FaceColor', [1, 1, 1] * 0.8, 'EdgeColor', 'k');
        hold on;

        n = from:to;
        pmf = binopdf(n, N, p(i));

        h_hist = histogram(data(:,i), ...
                           'Normalization', 'pdf', ...
                           'FaceColor', 'w', ....
                           'FaceAlpha', 1, ...
                           'LineWidth', line_width...
                           );

        plot(n, pmf, 'o', ...
             'MarkerSize', marker_size, ...
             'MarkerFaceColor', 'k', ...
             'MarkerEdgeColor', 'none');
        
        max_y = max([h_hist.Values pmf]) * 1.1;
        ylim([0 max_y]);
        xlim([from-0.5 to+0.5]);

        % Set tick labels and font sizes
        set(gca, ...
            'TickLabelInterpreter', 'latex', ...
            'FontSize', tickFontSize...
            );

        ax = gca;
        ax.YAxis.FontSize = tickFontSize;
        ax.XAxis.FontSize = tickFontSize;

        % X labels: Show parameter labels on the bottom row
        % X labels: Show parameter labels on the bottom row
        xlabel(labels{i}, ...
               'Interpreter', 'latex', ...
               'FontSize', legendFontSize ...
               );

        % Get the current x-tick values
        currentXTicks = xticks;
        
        if (max(currentXTicks) <= 5)
            xticks(0:max(currentXTicks));
        else
            % Retain only the integer ticks
            integerXTicks = currentXTicks(currentXTicks == floor(currentXTicks));
            xticks(integerXTicks);
        end
    end
end