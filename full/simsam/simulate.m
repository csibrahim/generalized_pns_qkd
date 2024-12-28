function [C, D0, D1, l, a, b, x, e] = simulate(N, thetaA, thetaB, thetaE, varargin)
    %simulate: Simulates detection events in a quantum key distribution (QKD) 
    %          system under a generalized Photon Number Splitting (PNS) attack.
    %
    % Inputs:
    %    N         - Total number of pulses to simulate
    %    thetaA    - Alice’s parameter set [lambdas, alpha, dAB]
    %    thetaB    - Bob’s parameter set [pa0, pa1, pc0, pc1, pd0, pd1, pe].
    %    thetaE    - Eve’s parameter set [dAE, pEB, k, Delta].
    %    varargin  - Optional name-value pairs:
    %                  'buffer'  - Batch size for simulation (default: 1e4)
    %                  'print'   - Display simulation progress (default: true)
    %
    % Outputs:
    %    C   - 1 x (8 * Nl) array containing counts for all detection events 
    %          across pulse intensities and basis configurations.
    %    D0  - Binary vector of detection events at detector D0
    %    D1  - Binary vector of detection events at detector D1
    %    l   - Intensity index for each pulse
    %    a   - Alice's basis choice (0 or 1) for each pulse
    %    b   - Bob's basis choice (0 or 1) for each pulse
    %    x   - Alice's bit choice (0 or 1) for each pulse
    %    e   - Eve's interception flag (1 if intercepted, 0 otherwise)
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Set up input parser with default values
    p = inputParser;
    addParameter(p, 'buffer', 1e4);
    addParameter(p, 'print', true);

    % Parse input arguments
    parse(p, varargin{:});

    % Assign parsed values to variables
    buffer = p.Results.buffer;
    print = p.Results.print;

    % Unpack system parameters
    [lambdas, alpha, dAB] = deal(thetaA{:});
    [pa0, pa1, pc0, pc1, pd0, pd1, pe] = deal(thetaB{:});
    [dAE, pEB, k, Delta] = deal(thetaE{:});

    Nl = length(lambdas);

    % Preallocate Arrays
    C = false(1, Nl * 8);  % Counts for each intensity and basis configuration
    D0 = false(1, N);      % Detector D0 events
    D1 = false(1, N);      % Detector D1 events
    a = false(1, N);       % Alice's basis choice
    b = false(1, N);       % Bob's basis choice
    x = false(1, N);       % Alice's bit choice
    e = false(1, N);       % Eve's interception flag
    l = zeros(1, N);       % Intensity index for each pulse

    if pa0 > 0 || pa1 > 0
        a0 = false(1, N);  % After-pulsing for D0
        a1 = false(1, N);  % After-pulsing for D1
    end

    % Define a Bernoulli sampler as a special case of binomial sampling
    bernrnd = @(p, varargin) binornd(1, p, varargin{:});

    % Simulate Pulses in Chunks
    num_chunks = ceil(N / buffer);

    for chunk = 1:num_chunks
        
        % Indices and size of the chunk
        first = (chunk - 1) * buffer + 1;
        last = min(chunk * buffer, N);
        ind = first:last;
        chunk_size = last - first + 1;

        % Randomize intensity, basis, and bit choices
        l(ind) = randi(Nl, 1, chunk_size);      % Intensity indices for this chunk
        a(ind) = bernrnd(0.5, 1, chunk_size);   % Alice's basis choices
        b(ind) = bernrnd(0.5, 1, chunk_size);   % Bob's basis choices
        x(ind) = bernrnd(0.5, 1, chunk_size);   % Alice's bit choices
        e(ind) = bernrnd(Delta, 1, chunk_size); % Determine if Eve intercepts the pulse

        intercepted  = e(ind);

        % Simulate photon counts and transmission probabilities
        
        ns = poissrnd(lambdas(l(ind))); % Photons emitted by Alice

        pAB = dist2prob(dAB, alpha);    % Alice to Bob
        pAE = dist2prob(dAE, alpha);    % Alice to Eve

        % Not intercepted
        ns(~intercepted) = binornd(ns(~intercepted), pAB);  % Photons arriving to Bob

        % Intercepted Pulses
        ns(intercepted) = binornd(ns(intercepted), pAE);  % Photons arriving to Eve
        ns(intercepted) = max(ns(intercepted) - k, 0);    % Photons intercepted by Eve
        ns(intercepted) = binornd(ns(intercepted), pEB);  % Photons arriving to Bob

        % Beam splitter probabilities
        p0 = p_bs(0, x(ind), a(ind), b(ind), pe);

        n0 = binornd(ns, p0);  % Photons to D0
        n1 = ns - n0;          % Photons to D1

        % Simulate detector clicks (dark counts and signal detections)
        d0 = bernrnd(pd0, 1, chunk_size);   % Dark counts at D0
        s0 = binornd(n0, pc0, 1, chunk_size);  % Photons detected at D0

        D0(ind) = d0 | (s0 >= 1); % Dark count or Detected at least 1 photon

        d1 = bernrnd(pd1, 1, chunk_size);   % Dark counts at D1
        s1 = binornd(n1, pc1, 1, chunk_size);  % Photons detected at D1

        D1(ind) = d1 | (s1 >= 1); % Dark count or Detected at least 1 photon

        % If after-pulsing is possible
        if pa0 > 0 || pa1 > 0

            % Go over the chunk
            for i = first:last
                
                if i > 1 % Ignore the first pulse

                    % If previous pulse is a click and not an after-pulse, 
                    % then consider this pulse to be an click

                    if D0(i - 1) && ~a0(i - 1)
                        a0(i) = rand < pa0;
                    end

                    if D1(i - 1) && ~a1(i - 1)
                        a1(i) = rand < pa1;
                    end
                end

                % True clicks or an after-pulse
                D0(i) = D0(i) | a0(i);
                D1(i) = D1(i) | a1(i);
            end

        end

        % Update detection counts
        
        C = C + countClicks(Nl, D0(ind), D1(ind), l(ind), a(ind), b(ind));

        % Display progress
        if print
            progress(chunk, num_chunks, 'Simulating clicks');
        end
    end
end

function C = countClicks(Nl, D0, D1, l, a, b)
    %countClicks: Counts detection events for all intensities and bases.
    %
    % Inputs:
    %     Nl  - Number of pulse intensities
    %     D0  - Binary detection events at D0
    %     D1  - Binary detection events at D1
    %     l   - Intensity level indices
    %     a   - Alice's basis choices
    %     b   - Bob's basis choices
    %
    % Outputs:
    %     C   - 1 x (8 * Nl) array of detection counts, categorized by:
    %           - Four detection types (00, 01, 10, 11)
    %           - Matching/non-matching bases
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    C = zeros(Nl, 4, 2);  % Detection counts for each intensity and basis match

    % Match condition: 1 = non-matching, 2 = matching
    match = (a == b) + 1;

    for i = 1:Nl
        for m = 1:2
            C(i, 1, m) = sum(~D1 & ~D0 & (l == i) & (match == m));  % No clicks (00)
            C(i, 2, m) = sum(~D1 &  D0 & (l == i) & (match == m));  % Click on D0 (01)
            C(i, 3, m) = sum( D1 & ~D0 & (l == i) & (match == m));  % Click on D1 (10)
            C(i, 4, m) = sum( D1 &  D0 & (l == i) & (match == m));  % Click on both (11)
        end
    end

    % Reshape to 1 x (8 * Nl) vector
    C = [reshape(C(:,:,1), 1, []) reshape(C(:,:,2), 1, [])];
end
