% experiment2_validate_hmm: Simulates detection events across various intensities 
%                           and compares observed distributions of clicks and 
%                           errors with theoretical distributions derived from
%                           both the i.i.d. and HMM model.
%
% Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
% Licensed under the MIT License (see LICENSE file for full details).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup                                                                  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng(1);                    % Set the random seed for reproducibility

restoredefaultpath; clear; % Restore default class path and clear workspace
addpath(genpath('../.'));  % Add all subfolders to the search path

% File path for loading/saving results (can be empty [])
file_path = 'saved_data/experiment2_validate_hmm'; 

% Set to 'true' to load data, 'false' to save after computation
loadData = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Session Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

runs = 1e4; % Total number of simulation runs

% Pulse options
N = 1e4;    % Number of pulses in each simulation run
Nl = 2;     % Number of intensity levels (lambdas)
limit = 10; % Maximum intensity

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bob's Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set detection parameters for Bob's system (based on GYS parameters)
pa = 0.1;     % After-pulsing probability
pc = 0.045;   % Detector efficiency
pd = 1.7e-6;  % Dark count probability
pe = 0.033;   % Misalignment probability

% Apply slight variation (+/- 10%) to represent detector differences
pa0 = pa * 0.9;   pa1 = pa * 1.1;
pc0 = pc * 0.9;   pc1 = pc * 1.1;
pd0 = pd * 0.9;   pd1 = pd * 1.1;

% Group Bob's parameters
thetaB = {pa0, pa1, pc0, pc1, pd0, pd1, pe}; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Alice's Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha = 0.21; % Attenuation coefficient in the channel
dAB = 50;     % Distance between Alice and Bob (in km)

% Get optimal intensity levels
lambdas = getLambdas(Nl, alpha, dAB, thetaB, limit);

% Group Alice's parameters
thetaA = {lambdas, alpha, dAB}; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Eve's Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dAE = 10;     % Distance between Alice and Eve
k = 3;        % Initial number of photons intercepted by Eve
Delta = 0.2;  % Initial interception rate by Eve

% Optimize pEB for Eve to remain undetected
pEB = get_pEB(thetaA, thetaB, dAE, k);

% Group Eve's parameters
thetaE = {dAE, pEB, k, Delta};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Priors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

max_k = getMaxK(thetaA, thetaB);
beta_k = 2 / (max_k - 1);


alphas = [1 1    1   2];
betas  = [2 1 beta_k 1];
ub = [dAB 1 inf 1];
lb = [0   0  1  0];

thetaP = {alphas,betas,ub,lb};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(~loadData)

    % Initialize result matrices for counts, errors, and detection events
    C = zeros(runs, 8 * Nl);     % Clicks per intensity/basis-match/detection-event 
    m = zeros(runs, 2 * Nl);     % Erroneous clicks per intensity/basis-match
    M = zeros(runs, 2 * Nl);     % XOR clicks per intensity/basis-match
    both = zeros(runs, 2 * Nl);  % AND clicks per intensity/basis-match
    
    % Run simulation across the specified number of runs
    for i = 1:runs
        
        % Simulate pulses and measure detection events
        [D0, D1, l, a, b, x] = simulate(N, thetaA, thetaB, thetaE, ...
                                        'print', false);

        % Count clicks for different basis match/non-match cases
        C(i,:) = countClicks(Nl, D0, D1, l, a, b);

        % Count signal and error clicks
        [m(i,:), M(i,:), both(i,:)] = measure(Nl, D0, D1, l, a, b, x);
    
        progress(i, runs, 'simulating sessions');
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load/Save Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if(~isempty(file_path))
        save(file_path, 'C', 'm', 'M', 'both');
    end

else
    load(file_path);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute HMM detection probabilities and plot results
Ps_hmm = P_hmm(thetaA, thetaB, thetaE);
plotPs(Ps_hmm, C);
set(gcf, 'Name', 'Probabilities - HMM', 'NumberTitle', 'off');

% Compute i.i.d. detection probabilities and plot results
Ps_iid = P_iid(thetaA, thetaB, thetaE);
plotPs(Ps_iid, C);
set(gcf, 'Name', 'Probabilities - i.i.d.', 'NumberTitle', 'off');

% Compute HMM error rates and gain, and visualize them
[deltaQs_hmm, Qs_hmm] = deltaQ_hmm(thetaA, thetaB, thetaE);
plot_deltaQ(N, m, M, deltaQs_hmm, Qs_hmm);
set(gcf, 'Name', 'Gain/Error - HMM', 'NumberTitle', 'off');

% Compute i.i.d. error rates and gain, and visualize them
[deltaQs_iid, Qs_iid] = deltaQ_iid(thetaA, thetaB, thetaE);
plot_deltaQ(N, m, M, deltaQs_iid, Qs_iid);
set(gcf, 'Name', 'Gain/Error - i.i.d.', 'NumberTitle', 'off');