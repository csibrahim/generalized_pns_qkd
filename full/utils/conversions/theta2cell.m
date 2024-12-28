function thetaCell = theta2cell(thetaR, varR)
    %theta2cell: Converts the matrix `thetaR` into a cell array `thetaCell`, while 
    %             handling the special case of the variable 'lambdas', which may span 
    %             multiple columns.
    %
    %             This function converts each column of `thetaR` into a cell entry, except
    %             for 'lambdas', which is extracted and stored as a single cell entry if
    %             present in `varR`.
    %
    % Inputs:
    %     thetaR    - Matrix of system parameters (n x m), where each column corresponds 
    %                 to a variable in `varR`
    %     varR      - Cell array of variable names corresponding to the columns of `thetaR`
    %
    % Outputs:
    %     thetaCell - Cell array where each entry corresponds to a variable in `varR`.
    %                 The 'lambdas' variable (if present) is combined into a single entry.
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Get the size of thetaR (n rows, m columns)
    [n, m] = size(thetaR);

    % Convert the matrix thetaR to a cell array where each column becomes a cell entry
    thetaCell = mat2cell(thetaR, n, ones(1, m));

    % Find the index of the variable 'lambdas' in varR
    [~, ind_lambdas] = ismember('lambdas', varR);

    % If 'lambdas' is present in varR
    if ind_lambdas
        % Compute the number of columns corresponding to 'lambdas'
        Nl = m - length(varR) + 1;

        % Extract the lambdas columns and combine them into a single cell entry
        lambdas = cell2mat(thetaCell(:, ind_lambdas:ind_lambdas + Nl - 1));

        % Remove the lambdas columns from thetaCell
        thetaCell(:, ind_lambdas:ind_lambdas + Nl - 1) = [];

        % Insert the combined lambdas entry back into thetaCell in the correct position
        thetaCell = [thetaCell(1:ind_lambdas-1), lambdas, thetaCell(ind_lambdas:end)];
    end

end