%experiment8_full_hmm: The fully Bayesian approach for simulating a session 
%                      (HMM) and inferring all parameters (none are fixed)
%
% Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
% Licensed under the MIT License (see LICENSE file for full details).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup                                                                  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Fix the random seed for reproducibility
seed = 1;

% File path for loading/saving results (can be empty [])
file_path = 'data/experiment8_full_hmm';

loadData = false;    % Set to 'true' to load data, 'false' to save after computation
savePlots = true;    % Set to 'true' to save plots

rng(seed);                 % Set the random seed
restoredefaultpath;        % Restore default path and clear workspace
addpath(genpath('../.'));  % Add all subfolders to the search path


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting Options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FigureWidth = 1200;  % Width of the figure in points
FontSize = 24;       % Font size for plot labels and text
CIp = 0.99;          % Confidence interval threshold

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Session Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define variable sets
variables = {'lambdas','alpha','dAB', ...
             'pa0','pa1','pc0','pc1','pd0','pd1','pe', ...
             'dAE','pEB','k','Delta'};

varR = variables; % Random variables (all)
varF = {};        % Fixed variables  (none)

% Variance options
noise = 0.05;      % Standard deviation of the noise to add during simulation
sigma = 2 * noise; % Standard deviation of the priors for random variables

% Pulse options
N = 1e9;    % Number of pulses in the simulation
Nl = 8;     % Number of intensity levels (lambdas)
limit = 10; % Maximum intensity

% Optimization options
algorithm = 'trust-region'; % Optimization algorithm, trust-region or quasi-newton
maxIters  = 1000;           % Maximum number of iterations

% Sampling options
Ns = 1e5;         % Number of samples
Nb = 1e3;         % Number of burn-in samples
display = false;  % Live histogram update (slow)
chunkSize = 10;   % Size of chunk per iteration
method = 'srss';  % Sampler options: 
                  %     'cmss'  - Covariance-Matching Slice Sampler
                  %     'srss'  - Shrinking-Rank Slice Sampler
                  %     'slice' - Slice sampling

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bob's Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Detection parameters for Bob's system (based on GYS parameters)
pa = 0.1;   % After-pulsing probability
pc = 0.045; % Detector efficiency
pd = 1.7e-6;% Dark count probability
pe = 0.033; % Misalignment probability

% Apply slight variation (+/- 10%) to represent detector differences
pa0 = pa * 0.9;   pa1 = pa * 1.1;
pc0 = pc * 0.9;   pc1 = pc * 1.1;
pd0 = pd * 0.9;   pd1 = pd * 1.1;

% Group Bob's parameters
thetaB = {pa0, pa1, pc0, pc1, pd0, pd1, pe}; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Alice's Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha = 0.21;     % Attenuation coefficient
dAB = 50;         % Distance between Alice and Bob (in km)

% Get optimal intensity levels
lambdas = get_lambdas(Nl, alpha, dAB, thetaB, limit);

% Group Alice's parameters
thetaA = {lambdas, alpha, dAB}; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Eve's Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dAE = 10;        % Distance between Alice and Eve
k = 3;           % Initial number of photons intercepted by Eve
Delta = 0.2;     % Initial interception rate by Eve

% Optimize pEB for Eve to remain undetected
pEB = get_pEB(thetaA, thetaB, dAE, k);

% Group Eve's parameters
thetaE = {dAE, pEB, k, Delta};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Collect parameters for Alice, Bob, and Eve
true_thetas = {thetaA, thetaB, thetaE};

% Process random and fixed variables
[varR, varF] = processVars(true_thetas, varF, varR, sigma); 

% Add noise to the random variables
sample_thetas = sampleParameters(true_thetas, varR, noise);

params = {true_thetas, sample_thetas, ...
          varF, varR, ...
          sigma, Nl, ...
          algorithm, maxIters, ...
          Ns, Nb, ...
          method, display, chunkSize};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (~loadData)
    % Simulate the session
    C = simulate(N, sample_thetas{:});

    % Assume random variables and HMM
    fprintf('\nAssuming random variables and HMM...\n');
    [thetaP, varR, thetaR, thetaR_MAP, samples, Ks, K_theory] = post_process(C, params);

    % Assume random variables and iid
    fprintf('\nAssuming random variables and i.i.d...\n');

    params_iid = params; % Copy original parameters

    % Set pa = 0 for true_theta
    params_iid{1}{2}{1} = 0;
    params_iid{1}{2}{2} = 0;
    
    [thetaP_iid, varR_iid, thetaR_iid, thetaR_MAP_iid, samples_iid, Ks_iid] = post_process(C, params_iid);

    % Assume fixed variables and HMM
    fprintf('\nAssuming fixed variables and HMM...\n');

    params_fixed = params; % Copy original parameters
    params_fixed{5} = 0;   % Set sigma = 0 to fix Alice's and Bob's parameters

    [thetaP_fixed, varR_fixed, thetaR_fixed, thetaR_MAP_fixed, samples_fixed, Ks_fixed] = post_process(C, params_fixed);

    % Assume fixed variables and iid
    fprintf('\nAssuming fixed variables and i.i.d...\n');

    params_fixed_iid = params_iid; % Copy iid parameters
    params_fixed_iid{5} = 0;       % Set sigma = 0 to fix Alice's and Bob's parameters

    [thetaP_fixed_iid, varR_fixed_iid, thetaR_fixed_iid, thetaR_MAP_fixed_iid, samples_fixed_iid, Ks_fixed_iid] = post_process(C, params_fixed_iid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save/Load Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    close all;
    
    if (~isempty(file_path))
        
        fprintf('Saving data ...\n');

        save(file_path, ...
             'samples', 'samples_fixed', 'samples_iid', 'samples_fixed_iid', ...
             'varR', 'varR_fixed', 'varR_iid', 'varR_fixed_iid', ...
             'thetaR', 'thetaR_fixed', 'thetaR_iid', 'thetaR_fixed_iid', ...
             'thetaP', 'thetaP_fixed',  'thetaP_iid',  'thetaP_fixed_iid', ...
             'thetaR_MAP', 'thetaR_MAP_fixed',  'thetaR_MAP_iid',  'thetaR_MAP_fixed_iid', ...
             'Ks', 'Ks_fixed','Ks_iid','Ks_fixed_iid', ...
             'K_theory');
    end

else
    
    fprintf('Loading data ...\n');
    load(file_path);

    fprintf('\nAssuming random variables and HMM...\n');
    printResults(thetaR_MAP, thetaR, varR);

    fprintf('\nAssuming random variables and i.i.d...\n');
    printResults(thetaR_MAP_iid, thetaR_iid, varR_iid);

    fprintf('\nAssuming fixed variables and HMM...\n');
    printResults(thetaR_MAP_fixed, thetaR_fixed, varR_fixed);

    fprintf('\nAssuming fixed variables and i.i.d...\n');
    printResults(thetaR_MAP_fixed_iid, thetaR_fixed_iid, varR_fixed_iid);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get LaTeX-formatted labels
[param_labels, K_labels] = getLabels(Nl, varR);
param_labels_fixed = getLabels(Nl, varR_fixed);
param_labels_iid = getLabels(Nl, varR_iid);
param_labels_fixed_iid = getLabels(Nl, varR_fixed_iid);

% Display posterior distributions of parameters assuming random variables and HMM
displayHistograms(samples, thetaP{3}, thetaP{4}, ...
                  'ground_truth', thetaR, ...
                  'labels', param_labels, ...
                  'CIp', CIp, ...
                  'FontSize', FontSize, ...
                  'FigureWidth', FigureWidth);
set(gcf, 'Name', 'Parameters - Full + HMM', 'NumberTitle', 'off');
full_hmm_params_plot = gcf;

% Display posterior distributions of parameters assuming fixed variables and HMM
displayHistograms(samples_fixed, thetaP_fixed{3}, thetaP_fixed{4}, ...
                  'ground_truth', thetaR_fixed, ...
                  'labels', param_labels_fixed, ...
                  'CIp', CIp, ...
                  'FontSize', FontSize, ...
                  'FigureWidth', FigureWidth);
set(gcf, 'Name', 'Parameters - Fixed + HMM', 'NumberTitle', 'off');
fixed_hmm_params_plot = gcf;

% Display posterior distributions of parameters assuming random variables and i.i.d.
displayHistograms(samples_iid, thetaP_iid{3}, thetaP_iid{4}, ...
                  'ground_truth', thetaR_iid, ...
                  'labels', param_labels_iid, ...
                  'CIp', CIp, ...
                  'FontSize', FontSize, ...
                  'FigureWidth', FigureWidth);
set(gcf, 'Name', 'Parameters - Full + i.i.d.', 'NumberTitle', 'off');
full_iid_params_plot = gcf;

% Display posterior distributions of parameters assuming fixed variables and i.i.d.
displayHistograms(samples_fixed_iid, thetaP_fixed_iid{3}, thetaP_fixed_iid{4}, ...
                  'ground_truth', thetaR_fixed_iid, ...
                  'labels', param_labels_fixed_iid, ...
                  'CIp', CIp, ...
                  'FontSize', FontSize, ...
                  'FigureWidth', FigureWidth);
set(gcf, 'Name', 'Parameters - Fixed + i.i.d.', 'NumberTitle', 'off');
fixed_iid_params_plot = gcf;

% Display distributions of secure key rates for each intensity assuming random variables and HMM
displayHistograms(Ks, ones(1, Nl), zeros(1, Nl), ...
                  'ground_truth', K_theory, ...
                  'labels', K_labels, ...
                  'CIp', CIp, ...
                  'FontSize', FontSize, ...
                  'FigureWidth', FigureWidth);
set(gcf, 'Name', 'Key Rates - Full + HMM', 'NumberTitle', 'off');
full_hmm_Ks_plot = gcf;

% Display distributions of secure key rates for each intensity assuming fixed variables and HMM
displayHistograms(Ks_fixed, ones(1, Nl), zeros(1, Nl), ...
                  'ground_truth', K_theory, ...
                  'labels', K_labels, ...
                  'CIp', CIp, ...
                  'FontSize', FontSize, ...
                  'FigureWidth', FigureWidth);
set(gcf, 'Name', 'Key Rates - Fixed + HMM', 'NumberTitle', 'off');
fixed_hmm_Ks_plot = gcf;

% Display distributions of secure key rates for each intensity assuming random variables and i.i.d.
displayHistograms(Ks_iid, ones(1, Nl), zeros(1, Nl), ...
                  'ground_truth', K_theory, ...
                  'labels', K_labels, ...
                  'CIp', CIp, ...
                  'FontSize', FontSize, ...
                  'FigureWidth', FigureWidth);
set(gcf, 'Name', 'Key Rates - Full + i.i.d.', 'NumberTitle', 'off');
full_iid_Ks_plot = gcf;

% Display distributions of secure key rates for each intensity assuming fixed variables and i.i.d.
displayHistograms(Ks_fixed_iid, ones(1, Nl), zeros(1, Nl), ...
                  'ground_truth', K_theory, ...
                  'labels', K_labels, ...
                  'CIp', CIp, ...
                  'FontSize', FontSize, ...
                  'FigureWidth', FigureWidth);
set(gcf, 'Name', 'Key Rates - Fixed + i.i.d.', 'NumberTitle', 'off');
fixed_iid_Ks_plot = gcf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save Figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if savePlots
    fprintf('Saving plots ...\n');
    
    print(full_hmm_params_plot, 'figures/full_hmm_params', '-dpdf');
    print(fixed_hmm_params_plot, 'figures/fixed_hmm_params', '-dpdf');
    print(full_iid_params_plot, 'figures/full_iid_params', '-dpdf');
    print(fixed_iid_params_plot, 'figures/fixed_iid_params', '-dpdf');

    print(full_hmm_Ks_plot, 'figures/full_hmm_Ks', '-dpdf');
    print(fixed_hmm_Ks_plot, 'figures/fixed_hmm_Ks', '-dpdf');
    print(full_iid_Ks_plot, 'figures/full_iid_Ks', '-dpdf');
    print(fixed_iid_Ks_plot, 'figures/fixed_iid_Ks', '-dpdf');

end