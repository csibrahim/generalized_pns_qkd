%experiment5_eve_iid: Simulate an iid session and infer the parameters of
%                     Eve, assuming all other parameters are fixed.
%
% Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
% Licensed under the MIT License (see LICENSE file for full details).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup                                                                  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Fix the random seed for reproducibility
seed = 1;

% File path for loading/saving results (can be empty [])
file_path = 'data/experiment5_eve_iid';

loadData = false;    % Set to 'true' to load data, 'false' to save after computation
savePlots = true;    % Set to 'true' to save plots

rng(seed);                 % Set the random seed
restoredefaultpath;        % Restore default path and clear workspace
addpath(genpath('../.'));  % Add all subfolders to the search path

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting Options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

scale = 2;                    % Factor scaling for screen readability
FigureWidth = 180;            % Width of the figure in mm
FontSize = 8;                 % Font size for plot labels and text
FontName = 'Times New Roman'; % Font name for plot labels and text
Interpreter = 'latex';        % Font rendering (latex or tex)
CIp = 0.99;                   % Confidence interval threshold

% Scale the figure size and font
FigureWidth = FigureWidth * scale;
FontSize = FontSize * scale;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Session Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
N = 1e9;    % Number of pulses in the simulation
Nl = 8;     % Number of intensity levels (lambdas)
limit = 10; % Maximum intensity

% Optimization options
algorithm = 'quasi-newton'; % Optimization algorithm, trust-region or quasi-newton
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

% Group parameters for Alice, Bob, and Eve
true_thetas = {thetaA, thetaB, thetaE};

% Process random and fixed variables
[varR, varF] = processVars(true_thetas, varF, varR, sigma); 

% Add noise to the random variables
sample_thetas = sampleParameters(true_thetas, varR, noise);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Group parameters
params = {true_thetas, sample_thetas, ...
          varF, varR, ...
          sigma, Nl, ...
          algorithm, maxIters, ...
          Ns, Nb, ...
          FontSize, FontName, Interpreter, FigureWidth, ...
          method, display, chunkSize};

if (~loadData)
    
    % Simulate the session
    C = simulate(N, true_thetas{:});

    % Sample and compute key rates
    [thetaP, varR, thetaR, thetaR_MAP, samples, Ks, K_theory] = post_process(C, params);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save/Load Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get labels
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save Figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if savePlots
    fprintf('Saving plots ...\n');
    print(params_plot, 'figures/eve_iid_params', '-dpdf');
    print(Ks_plot, 'figures/eve_iid_Ks', '-dpdf');
end