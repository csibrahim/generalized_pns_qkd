function [uppers, medians, lowers] = displayHistograms(samples, ub, lb, varargin)
    %displayHistograms: Generate histograms with kernel density estimates (KDE),
    %                  including confidence intervals and comparisons with ground truth values.
    %
    % Inputs:
    %   samples        - Matrix of samples (n_samples x n_variables)
    %   ub             - Upper bounds for each variable (1 x n_variables)
    %   lb             - Lower bounds for each variable (1 x n_variables)
    %   varargin       - Optional name-value pair arguments:
    %                      'ground_truth' - Ground truth values for each variable (default: [])
    %                      'labels'       - Cell array of labels for each variable (default: [])
    %                      'epsilon'      - Confidence level for computing CI (default: 1e-4)
    %                      'FontSize'     - Font size for labels and text (default: 30)
    %                      'FigureWidth'  - Width of the figure in points (default: 1035)
    %                      'newFigure'    - Create a new Figure (default: true)
    %
    % Outputs:
    %   uppers         - Upper bounds of the CI for each variable (1 x n_variables)
    %   medians        - Median values for each variable (1 x n_variables)
    %   lowers         - Lower bounds of the CI for each variable (1 x n_variables)
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Handle cell inputs for ub and lb if provided as cell arrays
    if iscell(ub), ub = [ub{:}]; end
    if iscell(lb), lb = [lb{:}]; end

    % Set up input parser with default values
    p = inputParser;
    addParameter(p, 'ground_truth', []);
    addParameter(p, 'labels', []);
    addParameter(p, 'epsilon', 1e-4);
    addParameter(p, 'FontSize', 30);
    addParameter(p, 'FigureWidth', 1035);
    addParameter(p, 'newFigure', true);

    % Parse input arguments
    parse(p, varargin{:});

    % Assign parsed values to variables
    FontSize = p.Results.FontSize;
    FigureWidth = p.Results.FigureWidth;
    ground_truth = p.Results.ground_truth;
    labels = p.Results.labels;
    epsilon = p.Results.epsilon;
    newFigure = p.Results.newFigure;

    if (iscell(ground_truth))
        ground_truth = [ground_truth{:}];
    end

    if newFigure
        prepFigure(FigureWidth, FontSize);
    end
    
    p = 1-epsilon;
    
    % Number of variables
    dim = length(ub);
    
    % Initialize vectors for CI bounds and medians
    uppers = zeros(1, dim);
    medians = zeros(1, dim);
    lowers = zeros(1, dim);

    % Number of samples
    n = size(samples, 1);

    % Legend labels
    label_mean = {'Mean'};
    label_KDE = {'KDE'};
    label_CI = {['$\mathrm{CI}_{', num2str(100 * p), '\%}$']};
    label_ground_truth = [];
    
    if ~isempty(ground_truth)
        label_ground_truth = 'Ground Truth';
    end
        
    
    % Set up subplots layout
    max_rows = ceil(sqrt(dim));
    cols = ceil(dim / max_rows);
    rows = min(max_rows, dim);

    
    labelFontSize = 1.5*FontSize / cols;
    legendFontSize = 0.75*FontSize;

    % Iterate through variables
    for i = 1:dim

        subplot(rows, cols, i);
        
        % Sort samples and compute empirical CDF
        x = sort(samples(:, i));
        cdf = (1:n) / n;
        
        % Compute confidence interval
        x_max = x(find(cdf >= 1 - epsilon, 1, 'first'));
        x_mid = x(find(cdf >= 0.5, 1, 'first'));
        x_min = x(find(cdf >= epsilon, 1, 'first'));
        
        if isempty(x_max)
            x_max = max(x);
        end

        if isempty(x_mid)
            x_mid = mean(x);
        end

        if isempty(x_min)
            x_min = min(x);
        end
        
        % Store CI bounds and medians
        uppers(i)  = x_max;
        medians(i) = x_mid;
        lowers(i)  = x_min;

        % Perform kernel density estimation
        if isinf(ub(i))
            
            x(x == lb(i)) = sqrt(realmin) + lb(i);
            
            [fs, xs] = ksdensity(x, ...
                                'Function', 'pdf', ...
                                'Support', 'positive', ...
                                'BoundaryCorrection', 'reflection'...
                                );
        else
            
            x(x == lb(i)) = sqrt(realmin) + lb(i);
            x(x == ub(i)) = (1 - eps) * (ub(i) - lb(i)) + lb(i);
            
            [fs, xs] = ksdensity(x, ...
                                'Function', 'pdf', ...
                                'Support', [lb(i), ub(i)], ...
                                'BoundaryCorrection', 'reflection'...
                                );
        end

        % Plot CI area
        top = max(fs) * 1.1;
        
        h_area = area([x_min, x_max], [1, 1] * top, ...
                      'FaceColor', [1, 1, 1] * 0.8, ...
                      'EdgeColor', 'k'...
                      );

        ys = [0, top];

        hold on;
        
        % Plot histogram
        histogram(x, ...
                  'Normalization', 'pdf', ...
                  'FaceColor', 'w', ...
                  'FaceAlpha', 1, ...
                  'EdgeColor', 'k', ...
                  'BinLimits', [x_min, x_max] ...
                  );

        % Plot KDE
        h_pdf = plot(xs, fs, 'k-', 'LineWidth', 1);

        % Plot ground truth if provided
        if ~isempty(ground_truth)
            xs = [ground_truth(i), ground_truth(i)];
            h_ground_truth = plot(xs, ys, 'r-', 'LineWidth', 2);
        else
            h_ground_truth = [];
        end

        % Plot mean
        mean_x = mean(x);
        xs = [mean_x, mean_x];
        h_mean = plot(xs, ys, 'k--', 'LineWidth', 2);

        ylim(ys);

        set(gca, 'TickLabelInterpreter', 'latex');
        
        % Display labels if provided
        if ~isempty(labels)
            xlabel(labels{i}, 'Interpreter', 'latex', 'FontSize', labelFontSize);
        end

        set(gca, 'YTick', [], 'YTickLabel', []);

        hold off;

    end

    % Store handlers for the legend
    legend_handlers = [h_ground_truth, h_mean, h_pdf, h_area];
    legend_labels = [label_ground_truth, label_mean, label_KDE, label_CI];

    % Adjust legend position to be horizontally centered
    subplot(rows, cols, 1);
    top_left = get(gca, 'Position');

    subplot(rows, cols, cols);
    top_right = get(gca, 'Position');

    top = top_left(2) + top_left(4);

    lh = legend(legend_handlers, legend_labels, 'Orientation', 'horizontal', 'Interpreter', 'latex', 'FontSize', legendFontSize);

    width = (top_right(1) + top_right(3) - top_left(1));
    center = top_left(1) + width / 2;

    x = center;
    y = (top + 1) / 2;
    newPosition = [x, y, 0, 0];

    set(lh, 'Position', newPosition);

end