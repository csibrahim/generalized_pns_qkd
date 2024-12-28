%opt_options: Defines default options for the fminunc function in
%             optimization routines used throughout the repository.
%
% Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
% Licensed under the MIT License (see LICENSE file for full details).

algorithm = 'quasi-newton';   % Quasi-Newton algorithm for optimization
tolX = 1e-30;                 % Termination tolerance for parameter values
tolFun = 1e-30;               % Termination tolerance for function value
optTol = 1e-30;               % Optimality tolerance for gradient
maxEval = 2e3;                % Max function evaluations allowed
maxIter = 1e2;                % Max iterations allowed

% Define optimization options for fminunc
options = optimoptions(...
    'fminunc', ...                         % Use fminunc solver
    'Algorithm', algorithm, ...            % Set optimization algorithm
    'TolX', tolX, ...                      % Parameter value tolerance
    'TolFun', tolFun, ...                  % Function value tolerance
    'OptimalityTolerance', optTol, ...     % Gradient tolerance
    'MaxFunctionEvaluations', maxEval, ... % Max function evaluations
    'MaxIterations', maxIter, ...          % Max iterations
    'Display', 'off', ...                  % Suppress display output
    'SpecifyObjectiveGradient', true);     % Gradients are provided
