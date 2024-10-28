function plotPs(Ps, C, varargin)
    % plotPs: Visualizes probability mass functions (PMFs) for detection configurations 
    %         across multiple intensities, accounting for matching and non-matching 
    %         basis settings and detector states (00, 01, 10, 11).
    %
    % Inputs:
    %   Ps          - Theoretical probabilities for each configuration across intensities
    %   C           - Observed counts for each configuration and intensity
    %   varargin    - Optional name-value pair arguments:
    %                   'epsilon'     - Threshold for confidence interval (default: 1e-4)
    %                   'FontSize'    - Font size for labels and text (default: 30)
    %                   'FigureWidth' - Width of the figure in points (default: 1035)
    %
    % Outputs:
    %   fig         - figure handle
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Set up input parser with default values
    p = inputParser;
    addParameter(p, 'epsilon', 1e-4);       % Default: 1e-4
    addParameter(p, 'FontSize', 30);        % Default: 30
    addParameter(p, 'FigureWidth', 1035);   % Default: 1035

    % Parse input arguments
    parse(p, varargin{:});

    % Assign parsed values to variables
    epsilon = p.Results.epsilon;
    FontSize = p.Results.FontSize;
    FigureWidth = p.Results.FigureWidth;
    
    % Font sizes for legend and layout
    legendFont = FontSize * 0.5;
    N = sum(C(1, :)); % Total trials, assuming same N across variables
    Nl = length(Ps) / 8; % Number of intensity levels

    % Define layout for subplots
    rows = 4;
    cols = 2 * Nl + 1;

    % Define indices for each detection configuration
    inds = cell(1, rows);
    for i = 1:rows
        inds{i} = [1:Nl, (4 * Nl + 1):5 * Nl] + (i - 1) * Nl;
    end

    % Prepare figure with specified width
    prepFigure(FigureWidth);

    % Row labels for detector states
    row_labels = {'00', '01', '10', '11'};

    % Loop through each row (detection configuration)
    for i = 1:rows

        ind = inds{i}; % Index for current row configuration
        not_matched = ind(1:Nl); % Indices for non-matching basis
        matched = ind((Nl + 1):2 * Nl); % Indices for matching basis
    
        % Adjust font size based on subplot layout
        subFontSize = FontSize * sqrt(2 / cols);
    
        % Plot PMF for non-matching and matching basis
        plotPMFs(N, Ps(not_matched), epsilon, C(:, not_matched), rows, i, cols, 1, subFontSize);
        plotPMFs(N, Ps(matched), epsilon, C(:, matched), rows, i, cols, Nl + 2, subFontSize);
    
        % Label for current row (detector configuration)
        th = text(1, 0.5, row_labels{i}, "Interpreter", "latex", 'Units', 'normalized', ...
                  'HorizontalAlignment', 'left', 'FontSize', FontSize);

        % Adjust position of row label
        set(th, 'Units', 'points');
        pos = get(th, 'Position');
        set(th, 'Position', [pos(1) + FontSize / 2, pos(2:end)]);
    end

    % Add labels for matching and non-matching basis across the subplots
    plotBasisLabels(rows, cols, FontSize);

    % Add legends to the PMF plots
    plotPMFLegends(epsilon, rows, cols, legendFont);
end