function plotPMFLegends(epsilon, rows, cols, FontSize)
    %plotPMFLegends: Adds a horizontal legend for confidence interval, data PMF, 
    %               and theoretical PMF to the top of a grid of subplots.
    %
    % Inputs:
    %   epsilon    - Epsilon value used to calculate confidence interval percentage
    %   rows       - Number of rows in the subplot grid
    %   cols       - Number of columns in the subplot grid
    %   FontSize   - Font size for the legend text
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Calculate confidence interval percentage for display in the legend
    p = 1 - epsilon;

    % Get positions of the left-most and right-most subplots
    subplot(rows, cols, 1);
    left_left = get(gca, "Position");
    subplot(rows, cols, cols);
    right_right = get(gca, "Position");

    % Create a horizontal legend for confidence interval, PMF data, and theory
    lh = legend(['$\mathrm{CI}_{', num2str(100 * p), '\%}$'], 'PMF (Data)', 'PMF (Theory)', ...
                'Orientation', 'horizontal', 'Interpreter', 'latex', 'FontSize', FontSize, 'LineWidth', FontSize / 100);

    set(lh, 'Units', 'normalized');  % Set units to normalized for better positioning
    pos = get(lh, "Position");       % Get legend dimensions

    % Center the legend between the left-most and right-most subplots
    x = left_left(1) + (right_right(1) + right_right(3) - left_left(1)) / 2;
    width = pos(3);
    height = pos(4);
    x = x - width / 2;
    y = 1 - height * 1.25;            % Place legend slightly below the top edge
    newPosition = [x, y, width, height];

    set(lh, 'Position', newPosition); % Update legend position
end