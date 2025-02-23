function fig = prepFigure(FigureWidth, ratio, FontSize, FontName, Interpreter, x_label, y_label)
    %prepFigure: Configures a figure with specified dimensions, font size, 
    %            and optional axis labels, centering it on the screen.
    %
    % Inputs:
    %   FigureWidth  - Width of the figure in mm
    %   ratio        - Height-to-width ratio (default: 1)
    %   FontSize     - Font size for labels and ticks (default: 8)
    %   FontName     - Name of the font for labels and text (default: Times New Roman)
    %   Interpreter  - Font rendering, 'latex' or 'tex' (default: 'latex')
    %   x_label      - X-axis label (default: none)
    %   y_label      - Y-axis label (default: none)
    %
    % Outputs:
    %   fig - Handle to the created figure
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).


    % Set default values if not provided
    if (nargin < 2), ratio = 1; end
    if (nargin < 3), FontSize = 8; end
    if (nargin < 4), FontName = 'Times New Roman'; end
    if (nargin < 5), Interpreter = 'latex'; end
    if (nargin < 6), x_label = []; end
    if (nargin < 7), y_label = []; end

    tickFontSize = max(FontSize * (5 / 8), 5);

    % Get screen size for centering
    screenSize = get(0, 'ScreenSize');

    % Get screen pixel density
    dpi = get(0, 'ScreenPixelsPerInch');

    % Convert from mm to points
    FigureWidth = FigureWidth * (72 / 25.4);

    % Convert from points to pixels
    FigureWidth = FigureWidth * (dpi / 72);

    % Compute figure height based on width and aspect ratio
    FigureHeight = FigureWidth * ratio;
    
    % Create a new figure
    figure;
    fig = gcf;

    % Remove extra padding around the axes
    set(gca, 'LooseInset', [0, 0, 0, 0]);

    % Center the figure on the screen
    y = screenSize(4) / 2 - FigureHeight / 2;
    x = screenSize(1) + (screenSize(3) - FigureWidth) / 2;
    set(fig, 'Position', [x y FigureWidth FigureHeight]);

    % Configure paper size for exporting the figure
    set(fig, 'PaperUnits', 'points');
    set(fig, 'PaperSize', [FigureWidth FigureHeight]);
    set(fig, 'PaperPosition', [0 0 FigureWidth FigureHeight]);
    
    % Set tick label formatting and font size
    set(gca, ...
        'FontSize', tickFontSize, ...
        'FontName', FontName, ...
        'TickLabelInterpreter', Interpreter);

    % Add x-axis label if provided
    if ~isempty(x_label)

        xlabel(x_label, ...
               'FontSize', FontSize, ...
               'FontName', FontName, ...
               'Interpreter', Interpreter);
    end

    % Add y-axis label if provided, adjust alignment
    if ~isempty(y_label)
        
        if strcmp(Interpreter, 'latex')
            y_label = ['$',y_label,'$'];
        end

        ylabel(y_label, ...
               'FontSize', FontSize, ...
               'FontName', FontName, ...
               'Interpreter', Interpreter, ...
               'Rotation', 0, ...
               'HorizontalAlignment', 'right', ...
               'VerticalAlignment', 'middle');
    end

end