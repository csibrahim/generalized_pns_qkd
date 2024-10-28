function plotBasisLabels(rows, cols, FontSize)
    %plotBasisLabels: Adds 'a≠b' and 'a=b' labels below subplots to distinguish 
    %                 basis match cases in a grid of subplots, adjusting label sizes 
    %                 and positions for alignment.
    %
    % Inputs:
    %   rows     - Number of subplot rows in the figure
    %   cols     - Total number of columns in the subplot grid
    %   FontSize - Font size for the labels (in points)
    %
    % This function creates annotations and lines below the subplots for 'a≠b' and 'a=b'
    % basis match cases, aligned according to the layout of subplots.
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    Nl = (cols - 1) / 2;  % Compute the number of intensities based on total columns

    fig = gcf;
    fig.Units = 'points';      % Set figure units to points for accurate sizing
    figPos = fig.Position;
    figHeight = figPos(4);     % Height of the figure in points

    % Access specific subplots to determine layout positions in normalized units
    subplot(rows, cols, (rows - 1) * cols + 1);
    left_left = get(gca, 'Position');  % Position of the leftmost subplot

    subplot(rows, cols, rows * cols);
    right_right = get(gca, 'Position');  % Position of the rightmost subplot

    subplot(rows, cols, (rows - 1) * cols + Nl);
    left_right = get(gca, 'Position');  % Right position of the 'a≠b' region

    subplot(rows, cols, (rows - 1) * cols + Nl + 2);
    right_left = get(gca, 'Position');  % Left position of the 'a=b' region

    % Calculate widths for annotation boxes based on subplot positions
    width1_norm = (left_right(1) + left_right(3)) - left_left(1);  % Width for 'a≠b' label
    width2_norm = (right_right(1) + right_right(3)) - right_left(1);  % Width for 'a=b' label

    % Define x-positions for labels
    x1_norm = left_left(1);  % x-position for 'a≠b'
    x2_norm = right_left(1);  % x-position for 'a=b'

    % Convert FontSize to normalized units relative to figure height
    FontSize_norm = 1.5 * FontSize / figHeight;

    % Add annotation for 'a≠b' label below
    annotation('textbox', [x1_norm, 0, width1_norm, FontSize_norm], 'String', '$a\neq b$', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', FontSize, 'EdgeColor', 'none', ...
        'Interpreter', 'latex', 'Units', 'normalized');

    % Line above 'a≠b' label
    lineX1_norm = [x1_norm, x1_norm + width1_norm];
    lineY1_norm = [FontSize_norm, FontSize_norm]; % Position line at the bottom of the text
    annotation('line', lineX1_norm, lineY1_norm, 'LineWidth', 0.5, 'Color', [0.5, 0.5, 0.5]);

    % Add annotation for 'a=b' label below
    annotation('textbox', [x2_norm, 0, width2_norm, FontSize_norm], 'String', '$a=b$', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', FontSize, 'EdgeColor', 'none', ...
        'Interpreter', 'latex', 'Units', 'normalized');

    % Line below 'a=b' label
    lineX2_norm = [x2_norm, x2_norm + width2_norm];
    lineY2_norm = [FontSize_norm, FontSize_norm]; % Position line at the bottom of the text
    annotation('line', lineX2_norm, lineY2_norm, 'LineWidth', 0.5, 'Color', [0.5, 0.5, 0.5]);
end