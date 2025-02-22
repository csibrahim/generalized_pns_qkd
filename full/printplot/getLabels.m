function [param_labels, K_labels] = getLabels(Nl, varR, Interpreter)
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
    
    param_labels = {'\lambda','\alpha','d_{AB}',...
                    'p^{}_{a^{}_{0}}','p^{}_{a^{}_{1}}','p^{}_{c^{}_{0}}','p^{}_{c^{}_{1}}','p^{}_{d^{}_{0}}','p^{}_{d^{}_{1}}','p^{}_{e}',...
                    'd^{}_{AE}','p^{}_{EB}','k','\Delta'};
    
    % Find the indices in 'variables' that match with 'varR'
    [~, ind] = ismember(varR, variables);

    % Select the labels corresponding to the variables in 'varR'
    param_labels = param_labels(ind);

    % Handle the case where 'lambdas' is present in varR
    [~, ind_lambdas] = ismember('lambdas', varR);

    if ind_lambdas
        % Generate labels for each lambda value (e.g., λ_1, λ_2)
        lambda_is = arrayfun(@(i) ['\lambda^{}_{' num2str(i) '}'], 1:Nl, 'UniformOutput', false);

        % Insert the expanded lambda labels into the labels array
        param_labels = [param_labels(1:ind_lambdas-1), lambda_is, param_labels(ind_lambdas+1:end)];
    end

    K_labels = arrayfun(@(x) sprintf('K^{}_{%d}', x), 1:Nl, 'UniformOutput', false);

    if strcmp(Interpreter,'latex')
        param_labels = cellfun(@(s) ['$', s, '$'], param_labels, 'UniformOutput', false);
        K_labels = cellfun(@(s) ['$', s, '$'], K_labels, 'UniformOutput', false);
        
    end

end