function log_pdf = logDirichlet(x, alpha)
    
    if(isrow(x))
        x = x';
        alpha = alpha';
    end
    
    log_gamma_sum_alpha = gammaln(sum(alpha));
    sum_log_gamma_alpha = sum(gammaln(alpha));
    sum_alpha_log_x = sum((alpha - 1) .* log(x));
    
    log_pdf = log_gamma_sum_alpha - sum_log_gamma_alpha + sum_alpha_log_x;
end