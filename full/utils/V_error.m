function variance = V_error(N, Q, delta)
    %V_error: This function calculates the variance of the error rate using the Delta Method.
    %
    % Inputs:
    %     N        - Number of trials (e.g., number of pulses)
    %     Q        - Success probability for each trial
    %     delta    - Error probability, conditioned on success
    %
    % Outputs:
    %     variance - The computed variance of the error rate
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Calculate the expected number of successes (K)
    E_K = N * Q;
    
    % Calculate the expected number of errors given the number of successes (M | K)
    E_M_given_K = E_K .* delta;

    % Calculate the variance of the number of successes (Var(K))
    Var_K = N * Q .* (1 - Q);
    
    % Calculate the variance of the number of errors given the number of successes (Var(M | K))
    Var_M_given_K = E_K .* delta .* (1 - delta);

    % Calculate the variance of the error rate using the Delta Method
    variance = (Var_M_given_K ./ E_K.^2) + (Var_K .* E_M_given_K.^2 ./ E_K.^4);

end