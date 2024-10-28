function plot_deltaQ(N, m, M, deltaQs, Qs, varargin)
    %plot_deltaQ: Visualizes probability mass functions (PMFs) across multiple 
    %             intensity levels for both matching and non-matching basis 
    %             configurations, showing error and detection probabilities.
    %
    % Inputs:
    %   N           - Number of samples for PMF calculations
    %   m           - Observed counts of erroneous clicks; for each intensity, non-matching
    %                 and matching configurations are included if split = true 
    %                 (size = [samples, 2*Nl]) or separately if split = false 
    %                 (size = [samples, Nl])
    %   M           - Observed counts of all clicks; format matches m
    %   deltaQs     - Theoretical error probabilities for each intensity level
    %   Qs          - Theoretical detection probabilities for each intensity level
    %   varargin    - Optional name-value pair arguments:
    %                   'epsilon'     - Threshold for confidence interval (default: 1e-4)
    %                   'split'       - Flag for matching/non-matching basis (default: true)
    %                   'FontSize'    - Font size for labels and text (default: 30)
    %                   'FigureWidth' - Width of the figure in points (default: 1035)
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Set up input parser with default values
    p = inputParser;
    addParameter(p, 'epsilon', 1e-4);
    addParameter(p, 'split', true);
    addParameter(p, 'FontSize', 30);
    addParameter(p, 'FigureWidth', 1035);

    % Parse input arguments
    parse(p, varargin{:});

    % Assign parsed values to variables
    epsilon = p.Results.epsilon;
    split = p.Results.split;
    FontSize = p.Results.FontSize;
    FigureWidth = p.Results.FigureWidth;
  
    % Set font sizes and determine subplot layout based on intensity count
    legendFont = FontSize * 0.75;
    Nl = size(m, 2); % Total number of intensities, split across matching if necessary

    % Adjust columns and split Nl if matching/non-matching are to be separated
    cols = Nl;
    if split
        cols = Nl + 1;
        Nl = Nl / 2;
    end

    rows = 2;                 % Two rows: one for m, one for M
    data = {m, M};            % Observed data for number of errors (m) and number of clicks (M)
    Ps = {deltaQs, Qs};       % Theoretical probabilities for error (deltaQs) and detection (Qs)
    row_labels = {'m', 'M'};  % Labels for number of errors (m) and number of clicks (M)

    % Prepare figure with specified width
    prepFigure(FigureWidth);

    % Loop over the two cases (m and M)
    for i = 1:2
        
        % Adjust font size based on subplot layout
        subFontSize = FontSize * sqrt(2 / cols);
        
        % Plot PMF for the first Nl intensities
        plotPMFs(N, Ps{i}(1:Nl), epsilon, data{i}(:, 1:Nl), rows, i, cols, 1, subFontSize);

        % If split is true, plot the PMF for the next Nl intensities
        if split
            plotPMFs(N, Ps{i}(Nl+1:end), epsilon, data{i}(:, Nl+1:end), rows, i, cols, Nl + 2, subFontSize);
        end

        % Add row labels ('m' or 'M') on the right side of each row
        th = text(1, 0.5, row_labels{i}, "Interpreter", "latex", 'Units', 'normalized', ...
                  'HorizontalAlignment', 'left', 'FontSize', FontSize);
            
        % Adjust position of row labels for readability
        set(th, 'Units', 'points');
        pos = get(th, 'Position');
        set(th, 'Position', [pos(1) + FontSize / 2, pos(2:end)]);
        set(th, 'Units', 'normalized');
    end
    
    % If split, add labels for matching/non-matching basis columns
    if split
        plotBasisLabels(rows, cols, FontSize);
    end

    % Add legends for the PMF plots
    plotPMFLegends(epsilon, rows, cols, legendFont);
end