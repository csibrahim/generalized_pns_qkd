function [param_labels, K_labels] = getLabels(Nl, varR)
    %getLabels: Generates LaTeX-formatted labels for a given list 
    %            of variables `varR`. If 'lambdas' is present, it expands it into multiple 
    %            indexed labels based on the number of lambda values (Nl).
    %
    % Inputs:
    %     Nl    - Number of distinct lambda values (expansion size for 'lambdas')
    %     varR  - Cell array of variable names for which labels are required
    %
    % Outputs:
    %     param_labels - Cell array of LaTeX-formatted labels corresponding to `varR`
    %     R_labels     - Cell array of LaTeX-formatted labels corresponding to `R_lambda`
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).
    
    % Define the standard variable names and their corresponding LaTeX labels
    variables = {'lambdas','alpha','dAB',...
                 'pa0','pa1','pc0','pc1','pd0','pd1','pe',...
                 'dAE','pEB','k','Delta'};
    
    param_labels = {'λ','α','d_{AB}',...
                    'p_{a_0}','p_{a_1}','p_{c_0}','p_{c_1}','p_{d_0}','p_{d_1}','p_{e}',...
                    'd_{AE}','p_{EB}','k','Δ'};
    
    % Find the indices in 'variables' that match with 'varR'
    [~, ind] = ismember(varR, variables);

    % Select the labels corresponding to the variables in 'varR'
    param_labels = param_labels(ind);

    % Handle the case where 'lambdas' is present in varR
    [~, ind_lambdas] = ismember('lambdas', varR);

    if ind_lambdas
        % Generate labels for each lambda value (e.g., λ_1, λ_2)
        lambda_is = arrayfun(@(i) ['λ_{' num2str(i) '}'], 1:Nl, 'UniformOutput', false);

        % Insert the expanded lambda labels into the labels array
        param_labels = [param_labels(1:ind_lambdas-1), lambda_is, param_labels(ind_lambdas+1:end)];
    end

    K_labels = arrayfun(@(x) sprintf('K_{%d}', x), 1:Nl, 'UniformOutput', false);

end