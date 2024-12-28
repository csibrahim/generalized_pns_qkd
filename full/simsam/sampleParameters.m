function sample_thetas = sampleParameters(thetas, varR, noise)
    %simulate: Samples values for random variables in varR 
    %
    % Inputs:
    %     thetas  - Cell array of parameter sets: {thetaA, thetaB, thetaE}
    %     noise   - Level of noise to add to the random variables
    %
    % Outputs:
    %     sample_thetas  - Cell array of sampled parameter values after adding noise
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    
    % Unpack theta parameters
    [thetaA, thetaB, thetaE] = deal(thetas{:});
    [lambdas, alpha, dAB] = deal(thetaA{:});
    [pa0, pa1, pc0, pc1, pd0, pd1, pe] = deal(thetaB{:});

    % Sample values based on noise and whether they are random variables
    lambdas = sampleValues('lambdas', lambdas, noise, inf, 0, varR);
    alpha = sampleValues('alpha', alpha, noise, inf, 0, varR);
    dAB = sampleValues('dAB', dAB, noise, inf, 0, varR);

    thetaA = {lambdas, alpha, dAB};

    pa0 = sampleValues('pa0', pa0, noise, 1, 0, varR);
    pa1 = sampleValues('pa1', pa1, noise, 1, 0, varR);
    pc0 = sampleValues('pc0', pc0, noise, 1, 0, varR);
    pc1 = sampleValues('pc1', pc1, noise, 1, 0, varR);
    pd0 = sampleValues('pd0', pd0, noise, 1, 0, varR);
    pd1 = sampleValues('pd1', pd1, noise, 1, 0, varR);
    pe = sampleValues('pe', pe, noise, 1, 0, varR);

    thetaB = {pa0, pa1, pc0, pc1, pd0, pd1, pe};

    sample_thetas = {thetaA, thetaB, thetaE};

end

function samples = sampleValues(variable, mu, noise, ub, lb, varR)
    % Helper function to sample values for a given variable
    
    % If not a random variable or the noise is 0, then just return mu
    if (~ismember(variable, varR) || noise == 0)
        samples = mu;  % No noise added
    else
        if isinf(ub)
            % Gamma distribution for semi-bounded variables
            
            % sigma is a fraction of mu
            sigma = mu * noise;
            % shift mu
            mu = mu - lb;

            % Infer parameters from mu and sigma
            [alpha, beta] = gamma_parameters(mu, sigma.^2);

            % Sample from a Gamma distribution with parameters alpha and
            % beta, and shifted by lb
            samples = gamrnd(alpha, 1 ./ beta) + lb;
        else
            % Beta distribution for bounded variables
            width = ub - lb;

            % sigma is a fraction of mu
            sigma = mu * noise / width;

            % Shift by lb and scale by width
            mu = (mu - lb) / width;

            % Infer parameters from mu and sigma
            [alpha, beta] = beta_parameters(mu, sigma.^2);

            % Sample from a Beta distribution with parameters alpha and
            % beta, scaled by width and shifted by lb
            samples = betarnd(alpha, beta) * width + lb;
        end
    end
end