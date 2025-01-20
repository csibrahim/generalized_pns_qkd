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
    %                 'EQs_ref'     - Another theoretical error probabilities to compare with (default: none)
    %                 'Qs_ref'      - Another theoretical detection probabilities to compare with (default: none)
    %                 'LegendNames' - Names for the baseline and reference (default: 'Theory' if Ps_ref is empty, otherwise {'P','Q'} )
    %                 'CIp'         - Confidence interval threshold (default: 0.99)
    %                 'split'       - Flag to split bases configurations (default: true)
    %                 'FontSize'    - Font size for text and labels (default: 8)
    %                 'FontName'    - Name of the font for labels and text (default: Times New Roman)
    %                 'FigureWidth' - Width of figure in mm (default: 180)
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Set up input parser with default values
    p = inputParser;
    addParameter(p, 'EQs_ref', []);
    addParameter(p, 'Qs_ref', []);
    addParameter(p, 'LegendNames', []);
    addParameter(p, 'CIp', 0.99);
    addParameter(p, 'split', true);
    addParameter(p, 'FontSize', 8);
    addParameter(p, 'FontName', 'Times New Roman');
    addParameter(p, 'FigureWidth', 180);

    % Parse input arguments
    parse(p, varargin{:});

    % Assign parsed values to variables
    EQs_ref = p.Results.EQs_ref;
    Qs_ref = p.Results.Qs_ref;
    LegendNames = p.Results.LegendNames;
    CIp = p.Results.CIp;
    split = p.Results.split;
    FontSize = p.Results.FontSize;
    FontName = p.Results.FontName;
    FigureWidth = p.Results.FigureWidth;

    if(isempty(LegendNames))
        if(isempty(Qs_ref))
            LegendNames = {''};
        else
            LegendNames = {'P - ','Q - '};
        end
    else
        LegendNames = cellfun(@(x) [x ' - '], LegendNames, 'UniformOutput', false);
    end

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
    
    if(~isempty(Qs_ref))
        Ps_ref = {Qs_ref, EQs_ref}; % Reference theoretical detection (Qs) and error (EQs) probabilities
    end

    % Labels for rows (signal and error counts)
    row_labels = {'G', 'R'};

    % Layout dimensions
    blocks = 3;  % Tile width in blocks
    rows = 2;    % Two rows for signal and error data

    % Aditional row and column in case of a split
    rows = blocks * rows + split;
    cols = blocks * cols + split;


    % Prepare figure and layout
    prepFigure(FigureWidth, rows/cols, FontSize, FontName);
    
    tiledlayout(rows, cols, ...
                'TileSpacing', 'compact', ...
                'Padding', 'compact');

    % Generate x-axis labels for intensity levels
    x_labels = arrayfun(@(x) sprintf('λ_{%d}', x), 1:Nl, 'UniformOutput', false);
    
    % Loop through data types (signal (S) and error (R) counts)
    for i = 1:numel(data)
        
        % Plot PMFs for each intensity level
        for j=1:Nl
            nexttile([blocks, blocks]);
            if(isempty(Qs_ref))
                p = {Ps{i}(j)};
            else
                p = {Ps{i}(j), Ps_ref{i}(j)};
            end
            plotPMFs(N, p, CIp, data{i}(:, j), FontSize, FontName, x_labels{j});
        end

        
        % Handle matching/non-matching basis split
        if split

            % Add a gap tile
            nexttile([blocks, 1]);
            axis off;

            % Plot PMFs for matching basis
            for j=1:Nl
                nexttile([blocks, blocks]);

                if(isempty(Qs_ref))
                    p = {Ps{i}(Nl+j)};
                else
                    p = {Ps{i}(Nl+j), Ps_ref{i}(Nl+j)};
                end

                plotPMFs(N, p, CIp, data{i}(:, Nl+j), FontSize, FontName, x_labels{j});
            end

        end
        
        % Add row labels ('G' or 'R') on the right of each row
        th = text(1, 0.5, row_labels{i}, ...
                  'Units', 'normalized', ...
                  'HorizontalAlignment', 'left', ...
                  'FontSize', FontSize, ...
                  'FontName', FontName);

        % Adjust position of row labels by half FontSize
        set(th, 'Units', 'points');
        posistion = get(th, 'Position');
        newPosistion = [posistion(1) + FontSize / 2, posistion(2:end)];

        set(th, 'Position', newPosistion);
        set(th, 'Units', 'normalized');
    end

    % Add legends for CI, data, and theory
    plotPMFLegends(CIp, FontSize, FontName, LegendNames);
    
    % Add legends for basis match (a = b and a ≠ b), if split = true
    if split
        plotBasisLabels(blocks, FontSize, FontName);
    end

end