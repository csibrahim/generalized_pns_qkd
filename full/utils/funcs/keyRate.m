function K = keyRate(Delta, delta, Q)
    %keyRate: Computes the secure-key rate using the GLLP formula,
    %         accounting for Eve's interception (Delta), the error
    %         rate (delta), and the gain (Q).
    %
    % Inputs:
    %     Delta  - Proportion of Eve's intercepted pulses (tagged pulses)
    %     delta  - Error rate (probability of bit errors)
    %     Q      - Gain (probability of Bob detecting a signal)
    %
    % Outputs:
    %     K      - Secure key rate, accounting for errors and eavesdropping
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Compute secure key rate using the GLLP formula
    K = Q .* (1 - Delta - H2(delta) - (1 - Delta) .* H2(delta ./ (1 - Delta)));

    % Check validity of the secure key rate
    fail = (K < 0) | isnan(K) | ~isreal(K); 

    % Set key rate to 0 for invalid cases
    K(fail) = 0;

end
