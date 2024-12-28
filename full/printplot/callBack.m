function stop = callBack(optimValues, state, maxIters)
    %callBack: This function serves as a callback to monitor the progress
    %          of an optimization routine. It updates the progress display
    %          and checks whether to stop the optimization based on the 
    %          current iteration and state.
    %
    % Inputs:
    %     optimValues - Structure containing information about the optimization:
    %                     'iteration' - The current iteration number of the optimizer
    %     state       - Current state of the optimization, specified as a string:
    %                     'iter' - Called at each iteration
    %                     'done' - Called when the optimization is complete
    %     maxIters    - Maximum number of iterations for the optimization
    %
    % Outputs:
    %     stop        - Boolean flag indicating whether to stop the optimization
    %                   (default: false, meaning the optimization continues)
    
    % Initialize stop flag as false. This indicates that the optimization
    % should continue unless explicitly set to true.
    stop = false;
    
    % Extract the current iteration number from the optimization values.
    iter = optimValues.iteration;

    % Switch based on the current state of the optimization process.
    switch state
        % If the current state is 'iter', update progress for each iteration.
        case 'iter'
            % Call the progress function to update the progress display.
            % iter+1 is used to account for the 0-based indexing.
            progress(iter+1, maxIters+1, 'optimizing');
        
        % If the state is 'done', perform a final progress update.
        case 'done'
            % Check if the current iteration is less than the maximum allowed iterations.
            if iter < maxIters
                % Call the progress function to finalize the progress display.
                % This ensures that the progress is marked as complete.
                progress(maxIters+1, maxIters+1, 'optimizing');
            end
    end
end