function plotPMFLegends(CIp, FontSize, FontName, Interpreter, LegendNames)
    %plotPMFLegends: Adds a horizontal legend for confidence intervals, 
    %                data PMF, and theoretical PMF at the top of a tile grid.
    %
    % Inputs:
    %   CIp         - Confidence interval threshold (default: 0.99)
    %   FontSize    - Font size for legend text
    %   FontName    - Name of the font to be used
    %   Interpreter - Font rendering, 'latex' or 'tex'
    %   LegendNames - Text labels for the legends
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    
    % Define legend entries for CI, data PMF, and theoretical PMF
    if strcmp(Interpreter, 'latex')
        CI_string = cellfun(@(x) [x,'$\mathrm{CI}_{', num2str(100*CIp), '\%}$'], LegendNames, 'UniformOutput', false);
    else
        CI_string = cellfun(@(x) [x,'CI(', num2str(100*CIp), '%)'], LegendNames, 'UniformOutput', false);
    end
    data_string = 'Data Histogram';
    theory_string = cellfun(@(x) [x,'PMF'], LegendNames, 'UniformOutput', false);

    legend_string = [CI_string(:)', {data_string}, theory_string(:)'];
    
    % Create horizontal legend with LaTeX formatting
    lh = legend(legend_string, ...
                'Orientation', 'horizontal', ...
                'FontSize', FontSize, ...
                'LineWidth', FontSize / 100, ...
                'FontName', FontName, ...
                'Interpreter', Interpreter, ...
                'NumColumns', 3);

    % Place legend at the top of the tile layout
    lh.Layout.Tile = 'north';
end