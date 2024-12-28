function [R, S, B] = EQ_measure(Nl, D0, D1, l, a, b, x)
    %EQ_measure: Computes the total number of errors (R), detected signal 
    %            (S), and double-clicks (B). These are computed separately 
    %            for matching and non-matching basis choices between Alice 
    %            and Bob
    %
    % Inputs:
    %     Nl  - Number of intensities
    %     D0  - Binary vector of detection events at detector D0
    %     D1  - Binary vector of detection events at detector D1
    %     l   - Intensity index for each pulse
    %     a   - Alice's basis choice (0 or 1) for each pulse
    %     b   - Bob's basis choice (0 or 1) for each pulse
    %     x   - Alice's bit choice (0 or 1) for each pulse
    %
    % Outputs:
    %     R   - Number of errors for each intensity
    %     S   - Number of signal events (D0 XOR D1) for each intensity
    %     B   - Number of double-clicks ( D0 AND D1) for each intensity
    %
    % Copyright (c) 2024 Ibrahim Almosallam
    % Licensed under the MIT License (see LICENSE file for full details).

    % Identify matching and non-matching basis choices
    match = a == b;

    % Compute errors, signals, and simultaneous clicks for non-matching bases
    [R0, S0, B0] = measure(Nl, D0(~match), D1(~match), l(~match), x(~match));
    
    % Compute errors, signals, and simultaneous clicks for matching bases
    [R1, S1, B1] = measure(Nl, D0(match), D1(match), l(match), x(match));

    % Combine results from matching and non-matching bases
    R = [R0, R1];
    S = [S0, S1];
    B = [B0, B1];
end

function [R, S, B] = measure(Nl, D0, D1, l, x)
    %measure_m: Computes the total number of errors (R), detected signal 
    %           (S), and double-clicks (B).
    %
    % Inputs:
    %     Nl  - Number of intensities
    %     D0  - Binary vector of detection events at detector D0
    %     D1  - Binary vector of detection events at detector D1
    %     l   - Intensity index for each pulse
    %     a   - Alice's basis choice (0 or 1) for each pulse
    %     b   - Bob's basis choice (0 or 1) for each pulse
    %     x   - Alice's bit choice (0 or 1) for each pulse
    %
    % Outputs:
    %     R   - Number of errors for each intensity
    %     S   - Number of signal events (D0 XOR D1) for each intensity
    %     B   - Number of double-clicks ( D0 AND D1) for each intensity
    %
    % Copyright (c) 2024 Ibrahim Almosallam
    % Licensed under the MIT License (see LICENSE file for full details).

    % Calculate XOR (either D0 or D1 clicks) and AND (both D0 and D1 click)
    xor_clicks = xor(D0, D1);
    and_clicks = and(D0, D1);

    % Determine errors based on Alice's bit choice (x)
    errors = (x & D0) | (~x & D1);

    % Initialize output arrays
    R = zeros(1, Nl);
    S = zeros(1, Nl);
    B = zeros(1, Nl);

    % Loop through each intensity level and compute statistics
    for i = 1:Nl
        S(i) = sum(xor_clicks & (l == i));           % Number of signal events
        R(i) = sum(xor_clicks & (l == i) & errors);  % Number of errors
        B(i) = sum(and_clicks & (l == i));           % Number of simultaneous clicks
    end
    
end
