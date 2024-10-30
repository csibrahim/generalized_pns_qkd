%opt_options: This script defines default options for the fminunc function 
%              used in optimization routines throughout the repository
%
% Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
% Licensed under the MIT License (see LICENSE file for full details).

% Set optimization algorithm and tolerance parameters
algorithm = 'trust-region';   % Use the trust-region algorithm for optimization
tolX = 1e-30;                 % Termination tolerance on the parameter values
tolFun = 1e-30;               % Termination tolerance on the function value
optTol = 1e-30;               % Optimality tolerance for the gradient
maxEval = 2e3;                % Maximum number of function evaluations
maxIter = 1e2;                % Maximum number of iterations

% Define optimization options for fminunc
options = optimoptions(...
    'fminunc', ...                      % Use fminunc solver
    'Algorithm', algorithm, ...         % Set optimization algorithm to trust-region
    'TolX', tolX, ...                   % Set tolerance on the parameter values
    'TolFun', tolFun, ...               % Set tolerance on the function value
    'OptimalityTolerance', optTol, ...  % Set optimality tolerance (gradient tolerance)
    'MaxFunctionEvaluations', maxEval, ... % Maximum number of function evaluations
    'MaxIterations', maxIter, ...       % Maximum number of iterations
    'Display', 'off', ...               % Turn off display output
    'SpecifyObjectiveGradient', true);  % Indicate that gradients are provided