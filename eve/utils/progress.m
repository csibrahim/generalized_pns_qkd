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

    persistent str; % Persistent string for previous output
    persistent startTime; % Persistent start time for total elapsed
    persistent lastTimeStr; % Store the last time string for proper backspacing
    persistent lastPercentage; % Store the last percentage displayed
    persistent lastUpdateTime; % Store the time of the last update

    % Initialize/Reset at the beginning of each new process
    if index == 1 || isempty(startTime)
        str = [];
        startTime = tic;
        lastTimeStr = [];
        lastPercentage = -1; % Initialize to impossible value to force first update
        lastUpdateTime = 0; % Initialize last update time
    end

    % Calculate elapsed time
    elapsedSecs = toc(startTime);

    timeSinceLastUpdate = (elapsedSecs - lastUpdateTime);

    if timeSinceLastUpdate < 0.1 && lastPercentage >= 0 && index < total
        return;
    end

    hours = floor(elapsedSecs / 3600);
    minutes = floor(mod(elapsedSecs, 3600) / 60);
    seconds = floor(mod(elapsedSecs, 60));

    % Calculate estimated time remaining
    progress_ratio = index / total;
    if progress_ratio > 0
        totalEstimatedSecs = elapsedSecs / progress_ratio;
        remainingSecs = totalEstimatedSecs - elapsedSecs;
    else
        remainingSecs = 0;
    end

    % Convert remaining seconds to hours, minutes, seconds
    remHours = floor(remainingSecs / 3600);
    remMinutes = floor(mod(remainingSecs, 3600) / 60);
    remSeconds = floor(mod(remainingSecs, 60));

    % Format elapsed time string
    if hours > 0
        timeStr = sprintf('[Elapsed %02d:%02d:%02d', hours, minutes, seconds);
    elseif minutes > 0
        timeStr = sprintf('[Elapsed %02d:%02d', minutes, seconds);
    else
        timeStr = sprintf('[Elapsed %d s', seconds);
    end

    % Format remaining time string
    if remHours > 0
        timeStr = sprintf('%s, Remaining %02d:%02d:%02d]', timeStr, remHours, remMinutes, remSeconds);
    elseif remMinutes > 0
        timeStr = sprintf('%s, Remaining %02d:%02d]', timeStr, remMinutes, remSeconds);
    else
        timeStr = sprintf('%s, Remaining %d s]', timeStr, remSeconds);
    end

    % Calculate current percentage
    currentPercentage = (index / total) * 100;

    % Check if either percentage or time string has changed
    % AND at least 0.1 seconds have elapsed since the last update
    if (abs(currentPercentage - lastPercentage) >= 0.01 || elapsedSecs == floor(elapsedSecs)) || index == 1 || index == total
        if index > 1
            % Calculate number of backspaces: length of previous string + '%% ' (2 chars)
            numBackspaces = length(str) + 2 + length(lastTimeStr);
            fprintf(repmat('\b', 1, numBackspaces));
        end
        
        % Update the progress string
        str = sprintf('%s %.2f', msg, currentPercentage);
        fprintf([str, '%% ', timeStr]);
        
        % Update the persistent variables
        lastTimeStr = timeStr;
        lastPercentage = currentPercentage;
        lastUpdateTime = elapsedSecs; % Update the last update time
    end

    % Print a newline when the process is complete
    if index == total
        fprintf('\n');
    end
end