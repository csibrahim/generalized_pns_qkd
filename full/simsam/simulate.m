function [D0, D1, l, a, b, x, e, sample_thetas] = simulate(N, thetas, varR, varargin)
    %simulate: Simulates the detection events for quantum key distribution (QKD) 
    %          under a generalized PNS attack
    %
    % Inputs:
    %     N          - Total number of pulses to simulate
    %     thetas     - Cell array of parameter sets: {thetaA, thetaB, thetaE}
    %     varR       - Cell array of variable names for random parameters
    %     varargin   - Optional name-value pair arguments:
    %                   'noise'   - Level of noise to add to the random variables (default: 0)
    %                   'match'   - Boolean flag to force basis match (default: false)
    %                   'buffer'  - Buffer size for batch simulation (default: 1e4)
    %                   'print'   - Print progress (default: true)
    %
    % Outputs:
    %     D0             - Logical array indicating detection events at detector D0
    %     D1             - Logical array indicating detection events at detector D1
    %     l              - Intensity index for each pulse
    %     a              - Binary array representing Alice's basis choice (0 or 1) for each pulse
    %     b              - Binary array representing Bob's basis choice (0 or 1) for each pulse
    %     x              - Binary array representing Alice's bit choice (0 or 1) for each pulse
    %     e              - Binary array indicating whether Eve intercepted the pulse (1 for intercept, 0 for no intercept)
    %     sample_thetas  - Cell array of sampled parameter values after adding noise
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    %% Parse Optional Inputs Using inputParser
    p = inputParser;
    addParameter(p, 'noise', 0);
    addParameter(p, 'match', false);
    addParameter(p, 'buffer', 1e4);
    addParameter(p, 'print', true);

    % Parse input arguments
    parse(p, varargin{:});

    % Assign parsed values to variables
    noise = p.Results.noise;
    match = p.Results.match;
    buffer = p.Results.buffer;
    print = p.Results.print;

    % Unpack theta parameters
    [thetaA, thetaB, thetaE] = deal(thetas{:});
    [lambdas, alpha, dAB] = deal(thetaA{:});
    [pa0, pa1, pc0, pc1, pd0, pd1, pe] = deal(thetaB{:});
    [dAE, pEB, k, Delta] = deal(thetaE{:});

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

    Nl = length(lambdas);

    % Preallocate detection event arrays
    D0 = false(1, N);  % Detection events at detector D0
    D1 = false(1, N);  % Detection events at detector D1
     x = false(1, N);   % Alice's bit choice (0 or 1)
     a = false(1, N);   % Alice's basis choice (0 or 1)
     b = false(1, N);   % Bob's basis choice (0 or 1)
     e = false(1, N);   % Eve's interception flag (1 if intercepted, 0 otherwise)
     l = zeros(1, N);   % Intensity index for each pulse

    num_chunks = ceil(N / buffer);

    for chunk = 1:num_chunks
        
        first = (chunk - 1) * buffer + 1;
        last = min(chunk * buffer, N);
        chunk_size = last - first + 1;

        l(first:last) = randi(Nl, 1, chunk_size);  % Intensity indices for this chunk
        x(first:last) = randi([0 1], 1, chunk_size);  % Alice's bit choices
        a(first:last) = randi([0 1], 1, chunk_size);  % Alice's basis choices
        
        if match
            b(first:last) = a(first:last);  % Force basis match
        else
            b(first:last) = randi([0 1], 1, chunk_size);  % Random basis choice for Bob
        end

        eve = rand(1, chunk_size) < Delta;  % Determine if Eve intercepts the pulse
        e(first:last) = eve;

        % Simulate number of photons in the pulse (Poisson distributed)
        ns = poissrnd(lambdas(l(first:last)));

        % Compute probabilities of transmission
        pAB = dist2prob(dAB, alpha);  % Probability of transmission from Alice to Bob
        pAE = dist2prob(dAE, alpha);  % Probability of transmission from Alice to Eve

        % Adjust photon counts for eavesdropping and losses
        ns(~eve) = binornd(ns(~eve), pAB);  % Transmission to Bob
        ns(eve) = binornd(ns(eve), pAE);    % Transmission to Eve
        ns(eve) = max(ns(eve) - k, 0);      % Eve intercepts photons
        ns(eve) = binornd(ns(eve), pEB);    % Eve transmits to Bob

        % Compute probabilities for the beam splitter
        p0 = p_bs(0, x(first:last), a(first:last), b(first:last), pe);

        n0 = binornd(ns, p0);  % Photons in path 0
        n1 = ns - n0;          % Photons in path 1

        % Simulate dark counts and detector clicks for both detectors
        d0 = binornd(1,pd0);   % Dark counts at detector D0
        s0 = binornd(n0,pc0);  % Detection at D0

        D0(first:last) = d0 | (s0 > 0);  % Combine dark counts and signal detection

        d1 = binornd(1,pd1);  % Dark counts at detector D1
        s1 = binornd(n1,pc1); % Detection at D1
        
        D1(first:last) = d1 | (s1 > 0);  % Combine dark counts and signal detection

        % Display progress
        if print
            progress(chunk, num_chunks, 'simulating clicks');
        end
    end

    % Simulate after-pulsing if applicable
    if (pa0 > 0 || pa1 > 0)
        a0 = false(1, N);  % After-pulsing for D0
        a1 = false(1, N);  % After-pulsing for D1

        for i = 2:N

            % Consider after-pulsing, only if the previous is a click and
            % not an after-pulse
            if (D0(i - 1) && ~a0(i - 1))
                a0(i) = rand < pa0;
            end

            if (D1(i - 1) && ~a1(i - 1))
                a1(i) = rand < pa1;
            end

            % Either a click or an after-pulse
            D0(i) = D0(i) | a0(i);
            D1(i) = D1(i) | a1(i);

            % Display progress
            if print && (mod(i, buffer) == 0 || i == 2 || i == N)
                progress(i - 1, N - 1, 'simulating after-pulse');
            end
        end
    end

end

function samples = sampleValues(variable, mu, noise, ub, lb, varR)
    % Helper function to sample values for a given variable
    
    % If not a random variable or the noise is 0, then just return mu
    if (~ismember(variable, varR) || noise == 0)
        samples = mu;  % No noise added
    else
        if isinf(ub)
            % Gamma distribution for unbounded variables
            
            % sigma is a fraction of mu
            sigma = mu * noise;
            % shift mu
            mu = mu - lb;

            % Infer parameters from mu and sigma
            [alpha, beta] = gamma_parameters(mu, sigma);

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
            [alpha, beta] = beta_parameters(mu, sigma);

            % Sample from a Beta distribution with parameters alpha and
            % beta, scaled by width and shifted by lb
            samples = betarnd(alpha, beta) * width + lb;
        end
    end
end