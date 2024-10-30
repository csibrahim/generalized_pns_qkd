function Rs = R(Delta, delta, Q)
    %R: This function computes the secure key rate based on the GLLP (Gottesman-Lo, 
    %    Lütkenhaus-Preskill) formula. The key rate accounts for the fraction of 
    %    Eve's interception (Delta), the error rate (delta), and the gain (Q).
    %
    % Inputs:
    %     Delta  - Proportion of Eve's intercepted pulses (tagged pulses)
    %     delta  - Error rate (probability of detecting a bit error)
    %     Q      - Gain (probability that Bob detects a signal)
    %
    % Outputs:
    %     Rs     - Secure key rate after accounting for errors and eavesdropping
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Compute the secure key rate using the GLLP formula
    Rs = Q .* (1 - Delta - H2(delta) - (1 - Delta) .* H2(delta ./ (1 - Delta)));

    % Set conditions for valid secure key rate calculation
    pass = delta ./ (1 - Delta) < 0.5;  % Condition that delta/(1-Delta) must be < 0.5
    neg = Rs < 0;                       % Check if any computed key rate is negative

    % Set Rs to 0 for invalid cases
    Rs(~pass | neg) = 0;

end