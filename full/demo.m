%demo_full: This script sets up and simulates a quantum key distribution 
%      session, allowing estimation of the parameters associated with an
%      eavesdropper using Bayesian inference and MCMC techniques. It saves 
%      and loads data, performs optimization, and computes secure key rates
%
% Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
% Licensed under the MIT License (see LICENSE file for full details).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup                                                                  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Fix the random seed for reproducibility
seed = 1;

% File path for loading/saving results (can be empty [])
file_path = 'demo_data';

loadData = false;    % Set to 'true' to load data, 'false' to save after computation
savePlots = true;    % Set to 'true' to save plots

rng(seed);                 % Set the random seed
restoredefaultpath;        % Restore default path and clear workspace
addpath(genpath('.'));     % Add all subfolders to the search path

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting Options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

scale = 2;                    % Factor scaling for screen readability
FigureWidth = 180;              % Width of the figure in mm
FontSize = 8;                   % Font size for plot labels and text
FontName = 'Times New Roman';   % Font name for plot labels and text
Interpreter = 'latex';          % Font rendering (latex or tex)
CIp = 0.99;                     % Confidence interval threshold

% Scale the figure size and font
FigureWidth = FigureWidth * scale;
FontSize = FontSize * scale;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Session Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define variable sets
variables = {'lambdas','alpha','dAB', ...
             'pa0','pa1','pc0','pc1','pd0','pd1','pe', ...
             'dAE','pEB','k','Delta'};

varR = {'dAE', 'pEB', 'k', 'Delta'};           % Random variables
varF = setdiff(variables, varR, 'stable');     % Fixed variables

% Variance options
noise = 0;   % Standard deviation of the noise to add during simulation
sigma = 0;   % Standard deviation of the priors for random variables

% Pulse options
N = 1e7;    % Number of pulses in the simulation
Nl = 8;     % Number of intensity levels (lambdas)
limit = 10; % Maximum intensity

% Optimization options
algorithm = 'quasi-newton'; % Optimization algorithm, trust-region or quasi-newton
maxIters  = 1000;           % Maximum number of iterations

% Sampling options
Ns = 1e4;         % Number of samples
Nb = 1e3;         % Number of burn-in samples
display = true;   % Live histogram update (slow)
chunkSize = 10;   % Size of chunk per iteration
method = 'slice'; % Sampler options: 
                  %     'cmss'  - Covariance-Matching Slice Sampler
                  %     'srss'  - Shrinking-Rank Slice Sampler
                  %     'slice' - Slice sampling

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bob's Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Detection parameters for Bob's system (based on GYS parameters)
pa = 0.0;   % After-pulsing probability
pc = 0.045; % Detector efficiency
pd = 1.7e-6;% Dark count probability
pe = 0.033; % Misalignment probability

% Apply slight variation (+/- 10%) to represent detector differences
pa0 = pa * 0.9;   pa1 = pa * 1.1;
pc0 = pc * 0.9;   pc1 = pc * 1.1;
pd0 = pd * 0.9;   pd1 = pd * 1.1;

% Group Bob's parameters
thetaB = {pa0, pa1, pc0, pc1, pd0, pd1, pe}; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Alice's Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha = 0.21;     % Attenuation coefficient
dAB = 50;         % Distance between Alice and Bob (in km)

% Get optimal intensity levels
lambdas = get_lambdas(Nl, alpha, dAB, thetaB, limit);

% Group Alice's parameters
thetaA = {lambdas, alpha, dAB}; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Eve's Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dAE = 10;        % Distance between Alice and Eve
k = 3;           % Initial number of photons intercepted by Eve
Delta = 0.2;     % Initial interception rate by Eve

% Optimize pEB for Eve to remain undetected
pEB = get_pEB(thetaA, thetaB, dAE, k);

% Group Eve's parameters
thetaE = {dAE, pEB, k, Delta};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Group parameters for Alice, Bob, and Eve
true_thetas = {thetaA, thetaB, thetaE};

% Process random and fixed variables
[varR, varF] = processVars(true_thetas, varF, varR, sigma); 

% Add noise to the random variables
sample_thetas = sampleParameters(true_thetas, varR, noise);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (~loadData)
    
    % Simulate the session
    C = simulate(N, true_thetas{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process and Get Prior Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Process fixed and random variables
    [varR, varF] = processVars(true_thetas, varF, varR, sigma);

    % Derive prior parameters and split variables into random and fixed
    [thetaP, thetaR, thetaF, varR, varF] = processParams(true_thetas, ...
                                                         sample_thetas, ...
                                                         varF, varR, sigma);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Maximum a posteriori probability (MAP) estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Find MAP estimation for random variables
    thetaR_MAP = MAP(C, varR, varF, thetaF, thetaP, ...
                     'algorithm', algorithm, ...
                     'maxIters', maxIters);

    % Display MAP results compared to ground truth
    printResults(thetaR_MAP, thetaR, varR);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical Integration through MCMC Sampling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Generate LaTeX-formatted labels
    param_labels = getLabels(Nl, varR, Interpreter);

    % Perform MCMC sampling for posterior distribution of random variables
    samples = sample(method, thetaR_MAP, Ns, Nb, ...
                     C, varR, varF, thetaF, thetaP, ...
                     'ground_truth', thetaR, ...
                     'display', display, ...
                     'FontSize', FontSize, ...
                     'FontName', FontName, ...
                     'Interpreter', Interpreter, ...
                     'FigureWidth', FigureWidth, ...
                     'labels', param_labels, ...
                     'chunkSize', chunkSize);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation of Secure Key Rate and Display Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Calculate secure key rate, gain, and QBER using true parameters
    [EQs, Qs, Delta] = EQ(thetaR, thetaF, varR, varF);   % Error and gain probabilities
    delta = EQs ./ Qs;                                   % QBER using true parameters
    K_theory = keyRate(Delta, delta, Qs);                % Secure key rate (theoretical)

    % Calculate secure key rate, gain, and QBER using sampled parameters
    thetaRs = theta2cell(samples, varR);                 % Convert samples to cell format
    [EQs, Qs, Deltas] = EQ(thetaRs, thetaF, varR, varF); % Error and gain probabilities
    deltas = EQs ./ Qs;                                  % QBER using sampled parameters
    Ks = keyRate(Deltas, deltas, Qs);                    % Secure key rate from samples
    
    % Discard non-matching results for final key rates
    K_theory = K_theory(Nl + 1:end);
    Ks = Ks(:, Nl + 1:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save/Load Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    close all; % close all open windows
    
    fprintf('Saving data ...\n');
    if (~isempty(file_path))
        save(file_path, 'samples', 'thetaR', 'thetaP', 'thetaR_MAP', 'K_theory', 'Ks');
    end

else
    fprintf('Loading data ...\n');
    load(file_path);
    printResults(thetaR_MAP, thetaR, varR);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get LaTeX-formatted labels
[param_labels, K_labels] = getLabels(Nl, varR, Interpreter);

% Display posterior distributions of parameters
[upper_params, medians_params, lowers_params] = displayHistograms(samples, thetaP{3}, thetaP{4}, ...
                                                                  'ground_truth', thetaR, ...
                                                                  'labels', param_labels, ...
                                                                  'CIp', CIp, ...
                                                                  'FontSize', FontSize, ...
                                                                  'FontName', FontName, ...
                                                                  'Interpreter', Interpreter, ...
                                                                  'FigureWidth', FigureWidth);
params_plot = gcf;

% Display distributions of secure key rates for each intensity
[upper_Ks, median_Ks, lower_Ks] = displayHistograms(Ks, ones(1, Nl), zeros(1, Nl), ...
                                                    'ground_truth', K_theory, ...
                                                    'labels', K_labels, ...
                                                    'CIp', CIp, ...
                                                    'FontSize', FontSize, ...
                                                    'FontName', FontName, ...
                                                    'Interpreter', Interpreter, ...
                                                    'FigureWidth', FigureWidth);

Ks_plot = gcf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save Figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if savePlots
    fprintf('Saving plots ...\n');
    print(params_plot, 'demo_params', '-dpdf');
    print(Ks_plot, 'demo_Ks', '-dpdf');
end