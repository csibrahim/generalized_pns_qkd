function plotPs(Ps, C, varargin)
    %plotPs: Visualizes PMFs for detection configurations across intensities, 
    %        distinguishing matching and non-matching bases for detector states.
    %
    % Inputs:
    %   Ps          - Theoretical probabilities for each configuration and intensity
    %   C           - Observed counts for each configuration and intensity
    %   varargin    - Optional name-value pairs:
    %                   'CIp'         - Confidence interval threshold (default: 0.99)
    %                   'FontSize'    - Font size for labels and text (default: 16)
    %                   'FigureWidth' - Figure width in points (default: 800)
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).


    % Set up input parser with default values
    p = inputParser;
    addParameter(p, 'CIp', 0.99);
    addParameter(p, 'FontSize', 16);
    addParameter(p, 'FigureWidth', 800);

    % Parse input arguments
    parse(p, varargin{:});

    % Assign parsed values to variables
    CIp = p.Results.CIp;
    FontSize = p.Results.FontSize;
    FigureWidth = p.Results.FigureWidth;

    % Calculate total trials and number of intensity levels
    N = sum(C(1, :));    % Total trials, assuming same N across variables
    Nl = length(Ps) / 8; % Number of intensity levels

    % Define layout for tiles
    blocks = 3;       % Tile width in blocks
    rows = 4;         % Rows for 4 detector configurations (00, 01, 10, 11)
    cols = 2 * Nl;    % Columns for matching and non-matching bases

    % Define indices for each detection configuration
    inds = cell(1, rows);
    for i = 1:rows
        inds{i} = [1:Nl, (4 * Nl + 1):5 * Nl] + (i - 1) * Nl;
    end

    % Determine layout ratio after adding the gap block
    ratio = (blocks*rows+1)/(blocks*cols+1);

    % Prepare figure with specified width and ratio
    prepFigure(FigureWidth, ratio);

    % Labels for rows and intensities
    row_labels = {'00', '01', '10', '11'};
    x_labels = arrayfun(@(x) sprintf('$\\lambda_{%d}$', x), 1:Nl, 'UniformOutput', false);

    % Create tiled layout for plots
    tiledlayout(blocks*rows+1, blocks*cols+1, ...
                'TileSpacing', 'compact', ...
                'Padding', 'compact');
    
    % Plot PMFs for each detection configuration (00, 01, 10, 11)
    for i = 1:numel(row_labels)

        ind = inds{i};                  % Index for current row configuration
        not_matched = ind(1:Nl);        % Indices for non-matching basis
        matched = ind((Nl + 1):2 * Nl); % Indices for matching basis
    
        % Plot PMF for non-matching
        for j=1:Nl
            nexttile([blocks,blocks]);
            plotPMFs(N, Ps(not_matched(j)), CIp, C(:, not_matched(j)), FontSize, x_labels{j});
        end

        % Add a gap tile
        nexttile([blocks,1]);
        axis 'off';
        
        % Plot PMF for matching
        for j=1:Nl
            nexttile([blocks,blocks]);
            plotPMFs(N, Ps(matched(j)), CIp, C(:, matched(j)), FontSize, x_labels{j});
        end
    
        % Add row labels for detector configurations
        th = text(1, 0.5, row_labels{i}, ...
                 'Interpreter', 'latex', ...
                 'Units', 'normalized', ...
                 'HorizontalAlignment', 'left', ...
                 'VerticalAlignment', 'middle', ...
                 'FontSize', FontSize);

        % Adjust position of row labels by half FontSize
        set(th, 'Units', 'points');
        position = get(th, 'Position');
        newPosition = [position(1) + FontSize / 2, position(2:end)];

        set(th, 'Position', newPosition);
        set(th, 'Units', 'normalized');

    end

    % Add legends for CI, data, and theory
    plotPMFLegends(CIp, FontSize);

    % Add legends for basis match (a = b and a ≠ b), if split = true
    plotBasisLabels(blocks, FontSize);

end