function makeSquares()
    %makeSquares Reshapes subplots into squares

    % Get all the subplot axes
    axesHandles = findall(gcf, 'Type', 'axes');
    
    % Loop through each subplot and set its position
    for i = 1:length(axesHandles)
        % Get the current position of the subplot
        pos = get(axesHandles(i), 'Position');
        
        % Determine the new size z
        z = min(pos(3), pos(4));
        
        % Set the new position, keeping the center the same
        newPos = [pos(1) + (pos(3) - z) / 2, pos(2) + (pos(4) - z) / 2, z, z];
        set(axesHandles(i), 'Position', newPos);
    end
end