function printResults(estimate, ground_truth)
    %printResults: This function prints the Maximum A Posteriori (MAP) estimation 
    %               results alongside the ground truth values for the random variables.
    %
    % Inputs:
    %     estimate      - The estimated values from the MAP procedure
    %     ground_truth  - The true values of the parameters
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Flatten the ground_truth cell array
    ground_truth = [ground_truth{:}];

    
    % Print table header
    fprintf("-------------------------------------------------------------\n");
    fprintf("Maximum a posteriori probability (MAP) estimation\n");
    fprintf("-------------------------------------------------------------\n\n");

    % Define the table's column names (variable names) and row labels (MAP estimate and Ground Truth)
    columns = {'dAE','$pEB','k','Delta'};
    rows = {'MAP Estimation', 'Ground Truth'};
    
    % Create the table with the MAP estimates and ground truth values
    T = array2table([estimate;ground_truth]);
    T.Properties.VariableNames = columns;
    T.Properties.RowNames = rows;
    
    % Display the table
    disp(T);
    fprintf("-------------------------------------------------------------\n");
    
end