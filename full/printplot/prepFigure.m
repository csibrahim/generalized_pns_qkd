function fig = prepFigure(FigureWidth, ratio, FontSize, x_label, y_label)
    %prepFigure: Configures a figure with specified dimensions, font size, 
    %            and optional axis labels, centering it on the screen.
    %
    % Inputs:
    %   FigureWidth  - Width of the figure in points
    %   ratio        - Height-to-width ratio (default: 1)
    %   FontSize     - Font size for labels and ticks (default: 16)
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
    if (nargin < 3), FontSize = 16; end
    if (nargin < 4), x_label = []; end
    if (nargin < 5), y_label = []; end

    % Get screen size for centering
    screenSize = get(0, 'ScreenSize');

    % Get screen pixel density
    dpi = get(0, 'ScreenPixelsPerInch');

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
        'TickLabelInterpreter', 'latex', ...
        'FontSize', 0.5 * FontSize);

    % Add x-axis label if provided
    if ~isempty(x_label)
        xlabel(x_label, ...
               'Interpreter', 'latex', ...
               'FontSize', FontSize);
    end

    % Add y-axis label if provided, adjust alignment
    if ~isempty(y_label)
        ylabel(y_label, ...
               'Interpreter', 'latex', ...
               'FontSize', FontSize, ...
               'Rotation', 0, ...
               'HorizontalAlignment', 'right', ...
               'VerticalAlignment', 'middle');
    end

end