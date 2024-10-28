function printResults(estimate, ground_truth, varR)
    %printResults: This function prints the Maximum A Posteriori (MAP) estimation 
    %               results alongside the ground truth values for the random variables.
    %
    % Inputs:
    %     estimate      - The estimated values from the MAP procedure
    %     ground_truth  - The true values of the random variables
    %     varR          - Cell array of random variable names
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Flatten the ground_truth cell array
    ground_truth = [ground_truth{:}];

    % Handle the case where 'lambdas' is a variable in varR
    [~, ind_lambdas] = ismember('lambdas', varR);
    if ind_lambdas
        % Compute the number of lambda values
        Nl = length(ground_truth) - length(varR) + 1;
        % Generate variable names for the lambda values (e.g., lambda_1, lambda_2, ...)
        lambda_is = arrayfun(@(i) ['lambda_' num2str(i)], 1:Nl, 'UniformOutput', false);
        % Replace 'lambdas' in varR with the lambda_i variables
        varR = [varR(1:ind_lambdas-1), lambda_is, varR(ind_lambdas+1:end)];
    end

    % Print table header
    fprintf("-------------------------------------------------------------\n");
    fprintf("Maximum a posteriori probability (MAP) estimation\n");
    fprintf("-------------------------------------------------------------\n\n");
    
    % Define the table's column names (variable names) and row labels (MAP estimate and Ground Truth)
    columns = varR;
    rows = {'MAP Estimation', 'Ground Truth'};
    
    % Create the table with the MAP estimates and ground truth values
    T = array2table([estimate; ground_truth]);
    T.Properties.VariableNames = columns;
    T.Properties.RowNames = rows;
    
    % Display the table
    disp(T);
    fprintf("-------------------------------------------------------------\n");
    
end