function plotBasisLabels(blocks, FontSize)
    %plotBasisLabels: Adds 'a≠b' and 'a=b' labels below tiles to indicate 
    %                 matching and non-matching bases in a subplot grid.
    %
    % Inputs:
    %   blocks    - Width of a tile in blocks
    %   FontSize  - Font size for the labels (points)
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Define labels for bases match cases
    bases_labels = {'$a\neq b$', '$a=b$'};
    
    for i=1:numel(bases_labels)
        
        % Create a tile spanning 2*blocks in width
        nexttile([1,2*blocks]);

        % Draw a line separator
        plot([0 1], [1 1], 'k-');
        
        % Add the label in the center of the tile
        text(0.5, 0.5, bases_labels{i}, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'top', ...
        'Interpreter', 'latex', ...
        'Units', 'normalized', ...
        'FontSize', FontSize);

        % Set axes limits to ensure label fits
        xlim([0 1]);
        ylim([0 1]);

        % Add a gap tile after the first label
        axis off;

        % Add a gap tile
        if(i==1)
            nexttile
            axis off;
        end
    end


end