function stop = callBack(optimValues, state, maxIters)
    %callBack: Monitors the progress of an optimization routine. Updates the 
    %          progress display and determines whether to stop the optimization 
    %          based on the current iteration and state.
    %
    % Inputs:
    %   optimValues - Structure containing optimization progress information:
    %                 - 'iteration': Current iteration number of the optimizer
    %   state       - Current state of the optimization process, specified as:
    %                 - 'iter': Called at each iteration
    %                 - 'done': Called when the optimization is complete
    %   maxIters    - Maximum number of iterations for the optimization
    %
    % Outputs:
    %   stop        - Boolean flag to indicate whether to stop the optimization:
    %                 - false: Continue optimization (default)
    %                 - true: Stop optimization
    %
    % Copyright (c) 2024 Ibrahim Almosallam
    % Licensed under the MIT License (see LICENSE file for full details).

    % Initialize the stop flag to false (continue optimization by default)
    stop = false;

    % Extract the current iteration number from the optimization values
    iter = optimValues.iteration;

    % Handle the current state of the optimization process
    switch state
        case 'iter'
            % Update progress display during each iteration
            % Use iter+1 to account for zero-based indexing
            progress(iter + 1, maxIters + 1, 'optimizing');

        case 'done'
            % Finalize progress display if optimization completes
            if iter < maxIters
                % Mark the progress as complete
                progress(maxIters + 1, maxIters + 1, 'optimizing');
            end
    end
end
