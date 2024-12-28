function plotPMFLegends(CIp, FontSize)
    %plotPMFLegends: Adds a horizontal legend for confidence intervals, 
    %                data PMF, and theoretical PMF at the top of a tile grid.
    %
    % Inputs:
    %   CIp      - Confidence interval threshold (default: 0.99)
    %   FontSize - Font size for legend text
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    
    % Define legend entries for CI, data PMF, and theoretical PMF
    CI_string = ['$\mathrm{CI}_{', num2str(100*CIp), '\%}$'];
    data_string = 'PMF (Data)';
    theory_string = 'PMF (Theory)';

    legend_string = {CI_string, data_string, theory_string};
    
    % Create horizontal legend with LaTeX formatting
    lh = legend(legend_string, ...
                'Orientation', 'horizontal', ...
                'Interpreter', 'latex', ...
                'FontSize', FontSize, ...
                'LineWidth', FontSize / 100);

    % Place legend at the top of the tile layout
    lh.Layout.Tile = 'north';
end