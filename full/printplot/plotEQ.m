function plotEQ(N, R, S, EQs, Qs, varargin)
    %plotEQ: Visualizes PMFs for error and detection probabilities across 
    %        intensity levels, distinguishing matching/non-matching bases.
    %
    % Inputs:
    %   N         - Number of samples for PMF calculations
    %   R         - Observed erroneous click counts; includes matching and 
    %               non-matching configurations if split=true (size: [samples, 2*Nl]) 
    %               or separately if split=false (size: [samples, Nl])
    %   S         - Observed signal click counts; format matches R
    %   EQs       - Theoretical error probabilities for each intensity
    %   Qs        - Theoretical detection probabilities for each intensity
    %   varargin  - Optional name-value pairs:
    %                 'CIp'         - Confidence interval threshold (default: 0.99)
    %                 'split'       - Flag to split bases configurations (default: true)
    %                 'FontSize'    - Font size for text and labels (default: 16)
    %                 'FigureWidth' - Width of figure in points (default: 800)
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Set up input parser with default values
    p = inputParser;
    addParameter(p, 'CIp', 0.99);
    addParameter(p, 'split', true);
    addParameter(p, 'FontSize', 16);
    addParameter(p, 'FigureWidth', 800);

    % Parse input arguments
    parse(p, varargin{:});

    % Assign parsed values to variables
    CIp = p.Results.CIp;
    split = p.Results.split;
    FontSize = p.Results.FontSize;
    FigureWidth = p.Results.FigureWidth;

    cols = size(R, 2);
    
    % Number of intensities, split across matching if necessary

    if split
        Nl = cols/2;
    else
        Nl = cols;
    end

    % Define observed data and theoretical probabilities
    data = {S, R};  % Observed signal (S) and error (R) counts
    Ps = {Qs, EQs}; % Theoretical detection (Qs) and error (EQs) probabilities

    % Labels for rows (signal and error counts)
    row_labels = {'$\mathbf{G}$', '$\mathbf{R}$'};

    % Layout dimensions
    blocks = 5;  % Tile width in blocks
    rows = 2;    % Two rows for signal and error data

    % Aditional row and column in case of a split
    rows = blocks * rows + split;
    cols = blocks * cols + split;


    % Prepare figure and layout
    prepFigure(FigureWidth, rows/cols);
    
    tiledlayout(rows, cols, ...
                'TileSpacing', 'compact', ...
                'Padding', 'compact');

    % Generate x-axis labels for intensity levels
    x_labels = arrayfun(@(x) sprintf('$\\lambda_{%d}$', x), 1:Nl, 'UniformOutput', false);
    
    % Loop through data types (signal (S) and error (R) counts)
    for i = 1:numel(data)
        
        % Plot PMFs for each intensity level
        for j=1:Nl
            nexttile([blocks, blocks]);
            plotPMFs(N, Ps{i}(j), CIp, data{i}(:, j), FontSize, x_labels{j});
        end

        
        % Handle matching/non-matching basis split
        if split

            % Add a gap tile
            nexttile([blocks, 1]);
            axis off;

            % Plot PMFs for matching basis
            for j=1:Nl
                nexttile([blocks, blocks]);
                plotPMFs(N, Ps{i}(Nl+j), CIp, data{i}(:, Nl+j), FontSize, x_labels{j});
            end

        end
        
        % Add row labels ('G' or 'R') on the right of each row
        th = text(1, 0.5, row_labels{i}, ...
                  'Interpreter', 'latex', ...
                  'Units', 'normalized', ...
                  'HorizontalAlignment', 'left', ...
                  'FontSize', FontSize);

        % Adjust position of row labels by half FontSize
        set(th, 'Units', 'points');
        posistion = get(th, 'Position');
        newPosistion = [posistion(1) + FontSize / 2, posistion(2:end)];

        set(th, 'Position', newPosistion);
        set(th, 'Units', 'normalized');
    end

    % Add legends for CI, data, and theory
    plotPMFLegends(CIp, FontSize);
    
    % Add legends for basis match (a = b and a ≠ b), if split = true
    if split
        plotBasisLabels(blocks, FontSize);
    end

end