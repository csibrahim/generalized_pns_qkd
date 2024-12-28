function progress(index, total, msg)
    %progress: Displays the progress of a process as a percentage, along with 
    %          elapsed and estimated remaining time, and a custom message.
    %
    % Inputs:
    %   index  - Current step or iteration number
    %   total  - Total number of steps or iterations
    %   msg    - Custom message to display with progress information
    %
    % Copyright (c) 2024 Ibrahim Almosallam
    % Licensed under the MIT License (see LICENSE file for full details).

    % Set defaults for update conditions
    minTime = 1e-1;    % Minimum time interval between updates (seconds)
    precision = 1e-2;  % Minimum percentage change required for updates
    
    % Persistent variables for tracking progress
    persistent lastUpdateTime lastPercent lastTimeDisplay lastPrintedLength startTime

    % Initialize persistent variables at the start of a new loop
    if index == 1 || isempty(startTime)
        startTime = tic;        % Record start time
        lastUpdateTime = 0;     % Time of the last update
        lastPercent = 0;        % Percentage of the last update
        lastTimeDisplay = 0;    % Last recorded elapsed time
        lastPrintedLength = 0;  % Length of the last printed message
    end

    currentTime = toc(startTime); % Elapsed time since start in seconds

    % Skip update if time interval is too short and not the final iteration
    if currentTime - lastUpdateTime < minTime && index ~= total
        return;
    end

    % Calculate percentage completion
    percentComplete = (index / total) * 100;
    
    % Determine if updates are necessary
    percentChanged = abs(percentComplete - lastPercent) >= precision;
    timeDisplayChanged = floor(currentTime) ~= lastTimeDisplay;

    if ~(percentChanged || timeDisplayChanged || index ==1 || index == total)
        return;% Skip update if no significant changes occurred
    end

    % Erase previous line using backspace characters
    eraseStr = repmat('\b', 1, lastPrintedLength);
    fprintf(eraseStr);

    % Estimate remaining time
    if index > 0 && index < total
        estimatedTotalTime = currentTime * total / index;
        remainingTime = estimatedTotalTime - currentTime;
    else
        remainingTime = 0;
    end

    % Format elapsed and remaining time
    elapsedTimeStr = formatTime(currentTime);

    % Construct the display string based on whether it's the last iteration
    if index == total
        % On the last iteration, only display elapsed time
        displayStr = sprintf('%s %.2f%% [Elapsed %s]', msg, percentComplete, elapsedTimeStr);
    else
        % Format remaining time
        remainingTimeStr = formatTime(remainingTime);
        % Display both elapsed and remaining time
        displayStr = sprintf('%s %.2f%% [Elapsed %s, Remaining %s]', msg, percentComplete, elapsedTimeStr, remainingTimeStr);
    end

    % Print the progress line
    fprintf('%s', displayStr);

    % Update the length of the last printed string
    lastPrintedLength = length(displayStr);

    % Update persistent variables
    lastUpdateTime = currentTime;
    lastPercent = percentComplete;
    lastTimeDisplay = floor(currentTime);

    % Finalize progress on the last iteration
    if index == total
        fprintf('\n');
    end

    % Nested function to format time
    function timeStr = formatTime(t)
        % Converts time in seconds to a formatted string (hh:mm:ss)
        if t < 10
            timeStr = sprintf('%d s', floor(t));
        elseif t < 60
            timeStr = sprintf('%02d s', floor(t));
        elseif t < 3600
            minutes = floor(t / 60);
            seconds = floor(mod(t, 60));
            timeStr = sprintf('%02d:%02d', minutes, seconds);
        else
            hours = floor(t / 3600);
            minutes = floor(mod(t, 3600) / 60);
            seconds = floor(mod(t, 60));
            timeStr = sprintf('%02d:%02d:%02d', hours, minutes, seconds);
        end
    end
end
