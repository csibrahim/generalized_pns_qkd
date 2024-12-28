function printResults(estimate, ground_truth, varR)
    %printResults: Displays MAP estimation results alongside ground truth 
    %              values in a formatted table. Handles multiple variables, 
    %              including lambda parameters.
    %
    % Inputs:
    %   estimate      - Estimated values from the MAP procedure
    %   ground_truth  - True values of the random variables
    %   varR          - Cell array of random variable names
    %
    % Copyright (c) 2024 Ibrahim Almosallam
    % Licensed under the MIT License.

    % Flatten ground_truth if provided as a cell array
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

    % Compute relative differences as percentages
    relative_diff = 100*(estimate - ground_truth) ./ ground_truth;

    % Prepare table data
    data = [ground_truth(:), estimate(:), relative_diff(:)];

    % Define row and column  labels
    rows = varR;
    columns = {'Ground Truth', 'MAP Estimation', 'Relative Difference'};

    % Create the table
    T = array2table(data, ...
                    'VariableNames', columns, ...
                    'RowNames', rows);

    % Convert relative differences to percentage strings
    T.(3) = strcat(num2str(T.(3), '%.2f'), '%');

    % Print table header
    fprintf("──────────────────────────────────────────────────────────────────────\n");
    fprintf("Maximum a Posteriori Probability (MAP) Estimation\n");
    fprintf("──────────────────────────────────────────────────────────────────────\n\n");
    
    % Display the table
    disp(T);
    
    % Print footer
    fprintf("──────────────────────────────────────────────────────────────────────\n");
end
