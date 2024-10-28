function fig = prepFigure(FigureWidth, FontSize, x_label, y_label)
    %prepFigure: Prepare and format a figure with specified width, font size, and axis labels.
    %            This function adjusts the figure size, font size, and axes labels, 
    %            while centering the figure on the screen.
    %
    % Inputs:
    %   FigureWidth  - Desired width of the figure in points
    %   FontSize     - Font size for axis labels and tick marks (default: 20)
    %   x_label      - X-axis label (default: none)
    %   y_label      - Y-axis label (default: none)
    %
    % Outputs:
    %   fig - figure handle
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Set default values if not provided
    if (nargin < 2), FontSize = 20; end
    if (nargin < 3), x_label = []; end
    if (nargin < 4), y_label = []; end

    % Create a new figure
    figure;
    fig = gcf;

    % Define tick font size relative to the provided FontSize
    tickFontSize = 0.5 * FontSize;

    % Get screen size for centering the figure
    screenSize = get(0, 'ScreenSize');

    % Set the figure units to points
    set(fig, 'Units', 'points');
    pos = get(fig, 'Position');

    % Calculate the height based on the specified width while maintaining aspect ratio
    FigureHeight = pos(4) * (FigureWidth / pos(3));

    % Center the figure on the screen
    y = screenSize(4) / 2 - FigureHeight / 2;
    x = screenSize(1) + (screenSize(3) - FigureWidth) / 2;
    set(fig, 'Position', [x y FigureWidth FigureHeight]);

    % Set paper size and position for exporting figures
    set(fig, 'PaperUnits', 'points');
    set(fig, 'PaperSize', [FigureWidth FigureHeight]);
    set(fig, 'PaperPosition', [0 0 FigureWidth FigureHeight]);

    % Set tick label interpreter to LaTeX and define tick font size
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize);

    % Add x-axis label if provided
    if ~isempty(x_label)
        xlabel(x_label, 'Interpreter', 'latex', 'FontSize', FontSize);
    end

    % Add y-axis label if provided and adjust its position
    if ~isempty(y_label)
        ylabelHandle = ylabel(y_label, 'Interpreter', 'latex', 'FontSize', FontSize, 'Rotation', 0);

        % Adjust the position of the y-axis label for better alignment
        set(ylabelHandle, 'Units', 'points');
        currentPosition = get(ylabelHandle, 'Position');
        set(ylabelHandle, 'Position', currentPosition - [FontSize, 0, 0]); % Shift label left by FontSize points
    end
end