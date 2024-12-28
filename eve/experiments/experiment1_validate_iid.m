%experiment1_validate_iid: Simulates detection events across various intensities 
%                          and compares observed distributions of clicks and 
%                          errors with theoretical distributions derived from
%                          the i.i.d. model.
%
% Copyright (c) 2024 Ibrahim Almosallam
% Licensed under the MIT License (see LICENSE file for full details).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup                                                                  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Fix the random seed for reproducibility
seed = 1;

% File path for loading/saving results (can be empty [])
file_path = 'data/experiment1_validate_iid'; 

loadData  = false;   % Set to 'true' to load data, 'false' to save after computation
savePlots = true;    % Set to 'true' to save plots

rng(seed);                 % Set the random seed
restoredefaultpath;        % Restore default path and clear workspace
addpath(genpath('../.'));  % Add all subfolders to the search path

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting Options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FigureWidth = 800;   % Width of the figure in points
FontSize = 16;       % Font size for plot labels and text
CIp = 0.99;          % Confidence interval threshold

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Session Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Pulse options
N  = 1e4;    % Number of pulses in each simulation run
Nr = 1e4;    % Total number of simulation runs
Nl = 2;      % Number of intensity levels (lambdas)
limit = 10;  % Maximum intensity

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bob's Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set detection parameters for Bob's system (based on GYS parameters)
pa = 0.0;    % After-pulsing probability
pc = 0.045;  % Detector efficiency
pd = 1.7e-6; % Dark count probability
pe = 0.033;  % Misalignment probability

% Apply slight variation (+/- 10%) to represent detector differences
pa0 = pa * 0.9;   pa1 = pa * 1.1;
pc0 = pc * 0.9;   pc1 = pc * 1.1;
pd0 = pd * 0.9;   pd1 = pd * 1.1;

% Group Bob's parameters
thetaB = {pa0, pa1, pc0, pc1, pd0, pd1, pe}; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Alice's Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha = 0.21; % Attenuation coefficient
dAB = 50;     % Distance between Alice and Bob (in km)

% Get optimal intensity levels
lambdas = get_lambdas(Nl, alpha, dAB, thetaB, limit);

% Group Alice's parameters
thetaA = {lambdas, alpha, dAB}; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Eve's Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dAE = 10;     % Distance between Alice and Eve
k = 3;        % Number of photons intercepted by Eve
Delta = 0.2;  % Interception rate by Eve

% Optimize pEB for Eve to remain undetected
pEB = get_pEB(thetaA, thetaB, dAE, k);

% Group Eve's parameters
thetaE = {dAE, pEB, k, Delta};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Priors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find the maximum number of photons Eve can intercept
k_max = get_k_max(thetaA, thetaB);

% Define the shape parameter for Eve's photon interception prior
beta_k = 2 / (k_max - 1);

alphas = [1 1    1   2]; % Prior alpha parameters
betas  = [2 1 beta_k 1]; % Prior beta parameters
ub = [dAB 1 inf 1];      % Upper bounds for the random variables
lb = [0   0  1  0];      % Lower bounds for the random variables

% Group prior parameters
thetaP = {alphas, betas, ub, lb};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(~loadData)

    % Initialize result matrices for counts, errors, and detection events
    C = zeros(Nr, 8 * Nl);  % Clicks per intensity/match/detection
    R = zeros(Nr, 2 * Nl);  % Erroneous clicks per intensity/match
    S = zeros(Nr, 2 * Nl);  % XOR clicks per intensity/match
    B = zeros(Nr, 2 * Nl);  % AND clicks per intensity/match
    
    % Run simulation across the specified number of runs
    for i = 1:Nr
        
        % Simulate pulses and measure detection events
        [C(i,:), D0, D1, l, a, b, x] = simulate(N, thetaA, thetaB, thetaE, ...
                                       'print', false);

        % Count signal and error clicks
        [R(i,:), S(i,:), B(i,:)] = EQ_measure(Nl, D0, D1, l, a, b, x);
    
        progress(i, Nr, 'simulating sessions');
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load/Save Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if(~isempty(file_path))
        fprintf('Saving data ...\n');
        save(file_path, 'C', 'R', 'S', 'B');
    end

else
    fprintf('Loading data ...\n');
    load(file_path);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute iid detection probabilities and plot results
Ps = P(thetaA, thetaB, thetaE);
plotPs(Ps, C, ...
         'FigureWidth', FigureWidth, ...
         'FontSize', FontSize, ...
         'CIp', CIp);

set(gcf, 'Name', 'Probabilities', 'NumberTitle', 'off');
Ps_plot = gcf;

% Compute iid error rates and gain, and plot results
[EQs, Qs] = EQ(thetaA, thetaB, thetaE);
plotEQ(N, R, S, EQs, Qs, ...
         'FigureWidth', FigureWidth, ...
         'FontSize', FontSize, ...
         'CIp', CIp);

set(gcf, 'Name', 'Gain/Error', 'NumberTitle', 'off');
EQs_plot = gcf;

% Compute and visualize iid error rates and gain for the decoy-state protocol
e0 = 1/2; % Random background noise
q = 1/2;  % QKD efficiency

% Normal channel efficiency
pAB = dist2prob(dAB, alpha);

% Consider double clicks as signals and errors
S = S+B;
R = R+B;

% Decoy-state formulations (only matching bases)
Qs_decoy = q * (pd + 1 - exp(-lambdas * pAB * pc))/Nl;
EQs_decoy = q * (e0 * pd + pe * (1 - exp(-lambdas * pAB * pc)))/Nl;

% Visualize the results
plotEQ(N, R(:, Nl+1:end), S(:, Nl+1:end), EQs_decoy, Qs_decoy, ...
       'split', false, ...
       'FigureWidth', FigureWidth/2, ...
       'FontSize', FontSize, ...
       'CIp', CIp);

set(gcf, 'Name', 'Gain/Error - Decoy', 'NumberTitle', 'off');
decoy_plot = gcf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save Figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if savePlots
    fprintf('Saving plots ...\n');
    print(Ps_plot, 'figures/validate_iid_Ps', '-dsvg');
    print(EQs_plot, 'figures/validate_iid_EQs', '-dsvg');
    print(decoy_plot, 'figures/validate_iid_decoy', '-dsvg');
end
