function progress(index, total, msg)
    %progress: Displays the progress of a process as a percentage and elapsed time with a custom message.
    %
    %            This function uses persistent variables to update the progress and only
    %            updates when the displayed information changes. It shows progress in 
    %            percentage form and elapsed time.
    %
    % Inputs:
    %     index  - Current step or iteration number
    %     total  - Total number of steps or iterations
    %     msg    - Custom message to display with the progress percentage
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Set default values
    minTime = 1e-1;    % Minimum time to pass before updating the progress
    precision = 1e-2;  % The percentage precision to be displayed
    
    persistent lastUpdateTime lastPercent lastTimeDisplay lastPrintedLength startTime

    % Check if this is the start of a new loop
    if index == 1 || isempty(startTime)
        % Initialize persistent variables
        startTime = tic; % Record start time
        lastUpdateTime = 0; % Time of last update
        lastPercent = 0;
        lastTimeDisplay = 0;
        lastPrintedLength = 0;
    end

    currentTime = toc(startTime); % Elapsed time since start in seconds

    % Allow update if it's the last iteration
    if currentTime - lastUpdateTime < minTime && index ~= total
        return;
    end

    % Compute percentage complete
    percentComplete = (index / total) * 100;
    
    % Check if percentage or time display has changed enough
    percentChanged = abs(percentComplete - lastPercent) >= precision;
    timeDisplayChanged = floor(currentTime) ~= lastTimeDisplay;

    if ~(percentChanged || timeDisplayChanged || index ==1 || index == total)
        % Neither percentage nor time display has changed enough
        return;
    end

    % Erase previous line using backspace characters
    eraseStr = repmat('\b', 1, lastPrintedLength);
    fprintf(eraseStr);

    % Compute estimated remaining time
    if index > 0 && index < total
        estimatedTotalTime = currentTime * total / index;
        remainingTime = estimatedTotalTime - currentTime;
    else
        remainingTime = 0;
    end

    % Format elapsed time
    elapsedTimeStr = formatTime(currentTime);

    % Build the display string based on whether it's the last iteration
    if index == total
        % On the last iteration, only display elapsed time
        displayStr = sprintf('%s %.2f%% [Elapsed %s]', msg, percentComplete, elapsedTimeStr);
    else
        % Format remaining time
        remainingTimeStr = formatTime(remainingTime);
        % Display both elapsed and remaining time
        displayStr = sprintf('%s %.2f%% [Elapsed %s, Remaining %s]', msg, percentComplete, elapsedTimeStr, remainingTimeStr);
    end

    % Print the new progress line
    fprintf('%s', displayStr);

    % Update the last printed length
    lastPrintedLength = length(displayStr);

    % Update persistent variables
    lastUpdateTime = currentTime;
    lastPercent = percentComplete;
    lastTimeDisplay = floor(currentTime);

    % If the loop is complete, print a new line and reset
    if index == total
        fprintf('\n');
    end

    function timeStr = formatTime(t)
        % Formats time according to specified rules
        % t is in seconds
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
