function plotPMFs(N, p, CIp, data, FontSize, FontName, label)
    %plotPMFs: Creates a bar plot comparing observed data with theoretical 
    %          binomial probabilities, including confidence intervals.
    %
    % Inputs:
    %   N         - Total trials for the binomial distribution
    %   p         - Success probability for each trial
    %   CIp       - Confidence level (e.g., 0.95 for 95% CI)
    %   data      - Observed data values
    %   FontSize  - Font size for labels (default: 8)
    %   FontName  - Name of the font for labels and text (default: Times New Roman)
    %   label     - X-axis label (default: none)
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).
   
    % Set default values for optional inputs
    if (nargin < 5),  FontSize = 8; end
    if (nargin < 6),  FontName = 'Times New Roman'; end
    if (nargin < 7),  label = []; end

    % Adjust font size for tick labels
    tickFontSize = max(FontSize * (5 / 8), 5);

    % Compute confidence interval bounds
    epsilon = (1 - CIp) / 2;
    
    % Calculate marker size and line width based on plot dimensions
    current_units = get(gca, 'Units');
    set(gca, 'Units', 'points');
    position = get(gca, "Position");

    marker_size = position(3) / 20;
    line_width = position(3)/300;

    set(gca, 'Units', current_units);

    k = numel(p);
    colors = hsv(k-1);
    colors = [0 0 0; colors];

    

    % Compute histogram counts
    [counts, edges] = histcounts(data, 'Normalization', 'pdf');
    
    % Compute bin centers for plotting
    binCenters = edges(1:end-1) + diff(edges) / 2;

    pmfs = cell(1, k);
    lowers = cell(1, k);
    uppers = cell(1, k);
    ns = cell(1, k);
    froms = zeros(1, k);
    tos = zeros(1, k);

    top = 0;
    
    for i=1:k

        % Define range for plot limits
        lowers{i} = binoinv(epsilon, N, p{i});
        uppers{i} = binoinv(1 - epsilon, N, p{i});
        
        froms(i) = min([data;lowers{i}]);
        tos(i) = max([data;uppers{i}]);

        % Compute theoretical PMF
        ns{i} = froms(i):tos(i);
        
        pmfs{i} = binopdf(ns{i}, N, p{i});
    
        % Determain plotting upper-limit
        top = max([counts(:); pmfs{i}(:); top]);
    end

    top = top * 1.1;

    hold on;
    
    for i = 1:k
        % Plot confidence interval region
        area([lowers{i}-0.5, uppers{i}+0.5], [1, 1] * top, ...
             'FaceColor', colors(i,:), ...
             'FaceAlpha', 0.2, ...
             'EdgeColor', 'k');
    end

    % Plot the histogram using bar
    h_hist = bar(binCenters, counts, ...
                 'FaceColor', 'w', ...
                 'FaceAlpha', 1, ...
                 'LineWidth', line_width...
                 );

    for i = 1:k
        % Plot theoretical PMF as markers
        scatter(ns{i}, pmfs{i}, 40, 'filled', ...
             'MarkerFaceColor', colors(i,:), ...
             'MarkerFaceAlpha', 0.5, ...
             'MarkerEdgeColor', 'none');
    end
    

    % Set axis limits
    ylim([0 top]);
    xlim([min(froms)-0.5 max(tos)+0.5]);

    % Format tick labels and adjust font sizes
    set(gca, ...
        'FontSize', tickFontSize, ...
        'FontName', FontName ...
        );

    ax = gca;
    ax.YAxis.FontSize = tickFontSize;
    ax.XAxis.FontSize = tickFontSize;

    % Add x-axis label
    xlabel(label, ...
           'Interpreter', 'tex', ...
           'FontSize', FontSize, ...
           'FontName', FontName ...
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