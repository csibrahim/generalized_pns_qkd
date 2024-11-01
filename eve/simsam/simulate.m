function [D0, D1, l, a, b, x, e] = simulate(N, thetaA, thetaB, thetaE, varargin)
    %simulate: Simulates the detection events for quantum key distribution (QKD) 
    %          under a generalized PNS attack
    %
    % Inputs:
    %    N         - Total number of pulses to simulate
    %    thetaA    - Cell array of system parameters specific to the channel 
    %                between Alice and Bob, containing [lambdas, alpha, dAB].
    %    thetaB    - Cell array of system parameters for Bob's detectors, 
    %                containing [pa0, pa1, pc0, pc1, pd0, pd1, pe].
    %    thetaE    - Cell array of system parameters related to Eve’s
    %                eavesdropping, containing [dAE, pEB, k, Delta].
    %    varargin  - Optional name-value pair arguments:
    %                  'match'   - Boolean flag to force basis match (default: false)
    %                  'buffer'  - Buffer size for batch simulation (default: 1e4)
    %                  'print'   - Print progress (default: true)
    %
    % Outputs:
    %     D0             - Logical array indicating detection events at detector D0
    %     D1             - Logical array indicating detection events at detector D1
    %     l              - Intensity index for each pulse
    %     a              - Binary array representing Alice's basis choice (0 or 1) for each pulse
    %     b              - Binary array representing Bob's basis choice (0 or 1) for each pulse
    %     x              - Binary array representing Alice's bit choice (0 or 1) for each pulse
    %     e              - Binary array indicating whether Eve intercepted the pulse (1 for intercept, 0 for no intercept)
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    %% Parse Optional Inputs Using inputParser
    p = inputParser;
    addParameter(p, 'match', false);
    addParameter(p, 'buffer', 1e4);
    addParameter(p, 'print', true);

    % Parse input arguments
    parse(p, varargin{:});

    % Assign parsed values to variables
    match = p.Results.match;
    buffer = p.Results.buffer;
    print = p.Results.print;

    % Unpack theta parameters
    [lambdas, alpha, dAB] = deal(thetaA{:});
    [pa0, pa1, pc0, pc1, pd0, pd1, pe] = deal(thetaB{:});
    [dAE, pEB, k, Delta] = deal(thetaE{:});

    
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

        l(first:last) = randi(Nl, 1, chunk_size);     % Intensity indices for this chunk
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