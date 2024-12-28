function printResults(estimate, ground_truth)
    %printResults: Displays MAP estimation results alongside ground truth 
    %              values in a formatted table with relative differences.
    %
    % Inputs:
    %   estimate      - Estimated values from the MAP procedure
    %   ground_truth  - True parameter values
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Flatten the ground_truth cell array
    ground_truth = [ground_truth{:}];

    % Compute relative differences as percentages
    relative_diff = 100*(estimate - ground_truth) ./ ground_truth;

    % Prepare table data
    data = [ground_truth(:), estimate(:), relative_diff(:)];

    % Define row and column  labels
    rows = {'dAE','pEB','k','Delta'};
    columns = {'Ground Truth', 'MAP Estimation', 'Relative Difference'};
    
    % Create the table
    T = array2table(data, ...
                    'VariableNames', columns, ...
                    'RowNames', rows);

    % Convert relative differences to percentage strings
    T.(3) = strcat(num2str(T.(3), '%.2f'), '%');
    
    % Print table header
    fprintf("──────────────────────────────────────────────────────────────────────\n");
    fprintf("Maximum a posteriori probability (MAP) estimation\n");
    fprintf("──────────────────────────────────────────────────────────────────────\n\n");
    
    % Display the table
    disp(T);

    % Print footer
    fprintf("──────────────────────────────────────────────────────────────────────\n");
    
end