function plotPs(Ps, C, varargin)
    %plotPs: Visualizes PMFs for detection configurations across intensities, 
    %        distinguishing matching and non-matching bases for detector states.
    %
    % Inputs:
    %   Ps          - Theoretical probabilities for each configuration and intensity
    %   C           - Observed counts for each configuration and intensity
    %   varargin    - Optional name-value pairs:
    %                   'Ps_ref'      - Another theoretical probabilities to compare with (default: none)
    %                   'LegendNames' - Names for the baseline and reference (default: 'Theory' if Ps_ref is empty, otherwise {'P','Q'} )
    %                   'CIp'         - Confidence interval threshold (default: 0.99)
    %                   'FontSize'    - Font size for labels and text (default: 8)
    %                   'FontName'    - Name of the font for labels and text (default: Times New Roman)
    %                   'FigureWidth' - Figure width in mm (default: 180)
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).


    % Set up input parser with default values
    p = inputParser;
    addParameter(p, 'Ps_ref', []);
    addParameter(p, 'LegendNames', []);
    addParameter(p, 'CIp', 0.99);
    addParameter(p, 'FontSize', 8);
    addParameter(p, 'FontName', 'Times New Roman');
    addParameter(p, 'FigureWidth', 180);

    % Parse input arguments
    parse(p, varargin{:});

    % Assign parsed values to variables
    Ps_ref = p.Results.Ps_ref;
    LegendNames = p.Results.LegendNames;
    CIp = p.Results.CIp;
    FontSize = p.Results.FontSize;
    FontName = p.Results.FontName;
    FigureWidth = p.Results.FigureWidth;

    if(isempty(LegendNames))
        if(isempty(Ps_ref))
            LegendNames = {'Theory'};
        else
            LegendNames = {'P','Q'};
        end
    end

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
    prepFigure(FigureWidth, ratio, FontSize, FontName);

    % Labels for rows and intensities
    row_labels = {'00', '01', '10', '11'};
    x_labels = arrayfun(@(x) sprintf('λ_{%d}', x), 1:Nl, 'UniformOutput', false);

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

            if(isempty(Ps_ref))
                p = {Ps(not_matched(j))};
            else
                p = {Ps(not_matched(j)), Ps_ref(not_matched(j))};
            end

            plotPMFs(N, p, CIp, C(:, not_matched(j)), FontSize, FontName, x_labels{j});

        end

        % Add a gap tile
        nexttile([blocks,1]);
        axis 'off';
        
        % Plot PMF for matching
        for j=1:Nl
            nexttile([blocks,blocks]);
            
            if(isempty(Ps_ref))
                p = {Ps(matched(j))};
            else
                p = {Ps(matched(j)), Ps_ref(matched(j))};
            end


            plotPMFs(N, p, CIp, C(:, matched(j)), FontSize, FontName, x_labels{j});
        end
    
        % Add row labels for detector configurations
        th = text(1, 0.5, row_labels{i}, ...
                 'Units', 'normalized', ...
                 'HorizontalAlignment', 'left', ...
                 'VerticalAlignment', 'middle', ...
                 'FontSize', FontSize, ...
                 'FontName', FontName);

        % Adjust position of row labels by half FontSize
        set(th, 'Units', 'points');
        position = get(th, 'Position');
        newPosition = [position(1) + FontSize / 2, position(2:end)];

        set(th, 'Position', newPosition);
        set(th, 'Units', 'normalized');

    end

    % Add legends for CI, data, and theory
    plotPMFLegends(CIp, FontSize, FontName, LegendNames);

    % Add legends for basis match (a = b and a ≠ b), if split = true
    plotBasisLabels(blocks, FontSize, FontName);

end