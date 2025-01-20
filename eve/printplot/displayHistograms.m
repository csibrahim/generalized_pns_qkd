function [uppers, medians, lowers] = displayHistograms(samples, ub, lb, varargin)
    %displayHistograms: Generate histograms with KDE, confidence intervals, 
    %                   and comparisons with ground truth values.
    %
    % Inputs:
    %   samples   - Sample matrix (n_samples x n_variables)
    %   ub        - Upper bounds for variables (1 x n_variables)
    %   lb        - Lower bounds for variables (1 x n_variables)
    %   varargin  - Optional name-value pairs:
    %                 'rows'         - Number of rows to display (default: auto)
    %                 'ground_truth' - Ground truth values (default: [])
    %                 'labels'       - Variable labels (default: [])
    %                 'CIp'          - Confidence level (default: 0.99)
    %                 'numBins'      - Number of histogram bins (default: 50)
    %                 'FontSize'     - Font size for labels (default: 8)
    %                 'FontName'     - Font size for labels (default: Times New Roman)
    %                 'FigureWidth'  - Figure width in mm (default: 180)
    %                 'newFigure'    - Flag to create a new figure (default: true)
    %                 'message'      - Message to display (default: none)
    %
    % Outputs:
    %   uppers   - CI upper bounds for each variable (1 x n_variables)
    %   medians  - Median values for each variable (1 x n_variables)
    %   lowers   - CI lower bounds for each variable (1 x n_variables)
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Handle cell inputs for ub and lb if provided as cell arrays
    if iscell(ub), ub = [ub{:}]; end
    if iscell(lb), lb = [lb{:}]; end

    % Number of variables
    dim = size(samples, 2);

    % Set up input parser with default values
    p = inputParser;
    addParameter(p, 'rows', ceil(sqrt(dim/4)));
    addParameter(p, 'ground_truth', []);
    addParameter(p, 'labels', []);
    addParameter(p, 'CIp', 0.99);
    addParameter(p, 'numBins', 50);
    addParameter(p, 'FontSize', 8);
    addParameter(p, 'FontName', 'Times New Roman');
    addParameter(p, 'FigureWidth', 180);
    addParameter(p, 'newFigure', true);
    addParameter(p, 'message', []);

    % Parse input arguments
    parse(p, varargin{:});

    % Assign parsed values to variables

    rows = p.Results.rows;
    ground_truth = p.Results.ground_truth;
    labels = p.Results.labels;
    CIp = p.Results.CIp;
    numBins = p.Results.numBins;
    FontSize = p.Results.FontSize;
    FontName = p.Results.FontName;
    FigureWidth = p.Results.FigureWidth;
    newFigure = p.Results.newFigure;
    message = p.Results.message;

    if (iscell(ground_truth))
        ground_truth = [ground_truth{:}];
    end

    % Adjust font size for tick labels
    tickFontSize = max((5 / 8) * FontSize, 5);

    % For CI calculation
    epsilon = (1 - CIp) / 2;
    
    % Initialize vectors for CI bounds and medians
    uppers = zeros(1, dim);
    medians = zeros(1, dim);
    lowers = zeros(1, dim);

    % Legend labels
    label_mean = {'Mean'};
    label_KDE = {'KDE'};
    label_CI = {['CI(', num2str(100 * CIp), '%)']};
    label_ground_truth = [];
    
    if ~isempty(ground_truth)
        label_ground_truth = 'Ground Truth';
    end
        
    
    % Set up subplots layout
    cols = ceil(dim / rows);
    ratio = rows/cols;

    if newFigure
        prepFigure(FigureWidth, ratio, FontSize);
    end

    blocks = 1 + (~isempty(message))*2; % Use 3 blocks in case of a message

    tiledlayout(blocks * rows + ~isempty(message), blocks * cols, ...
               'Padding', 'tight', ...
               'TileSpacing','compact');

    
    % Iterate through variables
    for i = 1:dim

        nexttile([blocks, blocks]);
        
        % Sort samples and compute empirical CDF
        x = sort(samples(:, i));


        x(x == lb(i)) = sqrt(realmin) + lb(i);

        % Determine the support of the KDE
        if isinf(ub(i))
            support = 'positive';
        else
            x(x == ub(i)) = (1 - eps) * (ub(i) - lb(i)) + lb(i);
            
            support = [lb(i) ub(i)];
        end

        % Perform kernel density estimation
        [pdfs, x_pdf] = ksdensity(x, ...
                                'Function', 'pdf', ...
                                'Support', support, ...
                                'BoundaryCorrection', 'log'...
                                );

        [cdfs, x_cdf] = ksdensity(x, ...
                                'Function', 'cdf', ...
                                'Support', support, ...
                                'BoundaryCorrection', 'log'...
                                );
        
        % Store CI bounds and medians
        uppers(i)  = x_cdf(find(cdfs >= 1 - epsilon, 1, 'first'));
        medians(i) = x_cdf(find(cdfs >= 0.5, 1, 'first'));
        lowers(i)  = x_cdf(find(cdfs >= epsilon, 1, 'first'));

        % Compute histogram counts
        [counts, edges] = histcounts(x, numBins, ...
                             'Normalization', 'pdf');
        
        % Compute bin centers for plotting
        binCenters = edges(1:end-1) + diff(edges) / 2;
        
        % Plot CI area
        top = max([counts(:); pdfs(:)]) * 1.1;

        h_area = area([lowers(i), uppers(i)], [1, 1] * top, ...
                      'FaceColor', [1, 1, 1] * 0.8, ...
                      'EdgeColor', 'k'...
                      );

        ys = [0, top];

        hold on;

        % Plot the histogram using bar
        bar(binCenters, counts, ...
            'FaceColor', 'w', ...
            'FaceAlpha', 1, ...
            'EdgeColor', 'k', ...
            'BarWidth', 1);

        % Plot KDE
        h_pdf = plot(x_pdf, pdfs, 'k-', 'LineWidth', 1);

        % Plot ground truth if provided
        if ~isempty(ground_truth)
            x_pdf = [ground_truth(i), ground_truth(i)];
            h_ground_truth = plot(x_pdf, ys, 'r-', 'LineWidth', 2);
        else
            h_ground_truth = [];
        end

        % Plot mean
        mean_x = mean(x);
        x_pdf = [mean_x, mean_x];
        h_mean = plot(x_pdf, ys, 'k--', 'LineWidth', 2);

        ylim(ys);

        set(gca, 'FontSize', tickFontSize, 'FontName', FontName);
        
        % Display labels if provided
        if ~isempty(labels)
            xlabel(labels{i}, ...
                   'FontSize', FontSize, ...
                   'FontName', FontName);
        end

        % Format Axes
        set(gca, ...
            'YTick', [], ...
            'YTickLabel', [] ...
            );

        axis tight;
        
        ax = gca;       
        ax.XLim(2) = max(ax.XLim(2), eps); 
        ax.YLim(2) = top;

        hold off;

    end

    % Display the mesaage on the last row if provided
    if(~isempty(message))
        nexttile([1, blocks*cols]);

        text(0.5, 1, message, ...
             'HorizontalAlignment', 'center', ...
             'VerticalAlignment', 'top', ...
             'Units', 'normalized', ...
             'FontSize', FontSize, ...
             'FontName', FontName);

         % Set axes length = 1
        xlim([0 1]);
        ylim([0 1]);

        % Turn off axes
        axis off;

    end

    % Store handlers for the legend
    legend_handlers = [h_area, h_ground_truth, h_mean, h_pdf];
    legend_labels = [label_CI, label_ground_truth, label_mean, label_KDE];


    % Plot the legends
    lh = legend(legend_handlers, legend_labels, ...
                'Orientation', 'horizontal', ...
                'FontSize', FontSize, ...
                'FontName', FontName);

    % Place them north
    lh.Layout.Tile = 'north';

end