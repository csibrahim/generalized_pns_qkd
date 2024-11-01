function [m, M, both] = measure(Nl, D0, D1, l, a, b, x)
    %measure: This function computes the total number of errors, detected events, 
    %          and simultaneous clicks for both detectors based on whether Alice 
    %          and Bob's basis choices match or not.
    %
    % Inputs:
    %     Nl  - Number of different intensity levels
    %     D0  - Binary array indicating whether detector 0 clicked
    %     D1  - Binary array indicating whether detector 1 clicked
    %     l   - Array of pulse intensities
    %     a   - Binary input, Alice's basis choice for each pulse (0 or 1)
    %     b   - Binary input, Bob's basis choice for each pulse (0 or 1)
    %     x   - Binary input, Alice's bit choice for each pulse (0 or 1)
    %
    % Outputs:
    %     m    - Total number of errors for each intensity level
    %     M    - Total number of detected events (XOR of D0 and D1 clicks) for each intensity level
    %     both - Total number of simultaneous clicks (AND of D0 and D1) for each intensity level
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Identify matching and non-matching basis choices
    match = a == b;

    % Compute the measures separately for matching and non-matching cases
    [m0, M0, both0] = measure_m(Nl, D0(~match), D1(~match), l(~match), x(~match));
    [m1, M1, both1] = measure_m(Nl, D0(match), D1(match), l(match), x(match));

    % Combine results for matching and non-matching cases
    M = [M0 M1];
    m = [m0 m1];
    both = [both0 both1];
    
end

function [m, M, both] = measure_m(Nl, D0, D1, l, x)
    % measure_m: This function computes the number of errors, detected events, and
    %            simultaneous clicks for non-matching or matching basis choices at each intensity level.
    %
    % Inputs:
    %     Nl    - Number of different intensity levels
    %     D0    - Binary array indicating whether detector 0 clicked
    %     D1    - Binary array indicating whether detector 1 clicked
    %     l     - Array of pulse intensities
    %     x     - Binary input, Alice's bit choice for each pulse (0 or 1)
    %
    % Outputs:
    %     m      - Total number of errors for each intensity level
    %     M      - Total number of detected events (XOR of D0 and D1 clicks) for each intensity level
    %     both   - Total number of simultaneous clicks (AND of D0 and D1) for each intensity level
    %
    % Copyright (c) 2024 Ibrahim Almosallam

    % Calculate XOR (either D0 or D1 clicks) and AND (both D0 and D1 click) of the detectors
    xor_clicks = xor(D0, D1);
    and_clicks = and(D0, D1);
    
    % Determine errors based on Alice's bit choice (x)
    errors = (x & D0) | (~x & D1);
    
    % Initialize output arrays
    m = zeros(1, Nl);
    M = zeros(1, Nl);
    both = zeros(1, Nl);
    
    % Loop through each intensity level and compute the statistics
    for i = 1:Nl
           M(i) = sum(xor_clicks & l == i);           % Total number of XOR clicks
           m(i) = sum(xor_clicks & l == i & errors);  % Total number of errors
        both(i) = sum(and_clicks & l == i);           % Total number of AND clicks
    end
end