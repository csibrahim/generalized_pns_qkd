function log_pmf = logMultinomial(x, p)
    % This function computes the log of the Multinomial PMF.
    % x: vector of counts (non-negative integers)
    % p: vector of probabilities (must sum to 1)
    
    % Convert x and p to column vectors if they are row vectors
    if isrow(x)
        x = x';
        p = p';
    end
    
    
    % Compute total number of trials
    n = sum(x);
    
    % Compute log of Multinomial PMF
    log_gamma_n_plus_1 = gammaln(n + 1);           % log(Gamma(n + 1)) = log(n!)
    sum_log_gamma_x_plus_1 = sum(gammaln(x + 1));  % sum(log(Gamma(x_i + 1))) = sum(log(x_i!))
    sum_x_log_p = sum(x .* log(p));                % sum(x_i * log(p_i))
    
    log_pmf = log_gamma_n_plus_1 - sum_log_gamma_x_plus_1 + sum_x_log_p;
end