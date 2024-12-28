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

loadData  = false;   % Set to 'true' to load data, 'false' to save after computation
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

% Set detection parameters for Bob's system (based on GYS parameters)
pa = 0.0;     % After-pulsing probability
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
k = 3;        % Initial number of photons intercepted by Eve
Delta = 0.2;  % Initial interception rate by Eve

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

% Group parameters
params = {thetaA, thetaB, thetaE, thetaP, ...
          algorithm, maxIters, ...
          Ns, Nb, ...
          method, display, chunkSize};

if(~loadData)

    % Simulate the session
    C = simulate(N, thetaA, thetaB, thetaE);

    % Sample and compute key rates
    [thetaE_MAP, samples, Ks, K_theory] = post_process(C, params);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save/Load Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    close all; % close all open windows
    
    if(~isempty(file_path))
        fprintf('Saving data ...\n');
        save(file_path, 'samples', 'thetaP', 'thetaE', 'thetaE_MAP', 'K_theory', 'Ks');
    end

else
    fprintf('Loading data ...\n');
    load(file_path);
    printResults(thetaE_MAP, thetaE);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get LaTeX-formatted labels
param_labels = {'$d_{\mathrm{AE}}$','$p_{\mathrm{EB}}$','$k$','$\Delta$'};
K_labels = arrayfun(@(x) sprintf('$K_{%d}$', x), 1:Nl, 'UniformOutput', false);

% Display posterior distributions of parameters
[upper_params, medians_params, lowers_params] = displayHistograms(samples, thetaP{3}, thetaP{4}, ...
                                                                  'ground_truth', thetaE, ...
                                                                  'labels', param_labels, ...
                                                                  'CIp', CIp, ... 
                                                                  'FontSize', FontSize, ...
                                                                  'FigureWidth', FigureWidth);
params_plot = gcf;

% Display distributions of secure key rates for each intensity
[upper_Ks, median_Ks, lower_Ks] = displayHistograms(Ks, ones(1,Nl), zeros(1,Nl), ...
                                                    'ground_truth', K_theory, ...
                                                    'labels', K_labels, ...
                                                    'CIp', CIp, ... 
                                                    'FontSize', FontSize, ...
                                                    'FigureWidth', FigureWidth);

Ks_plot = gcf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save Figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if savePlots
    fprintf('Saving plots ...\n');
    print(params_plot, 'figures/eve_iid_params', '-dsvg');
    print(Ks_plot, 'figures/eve_iid_Ks', '-dsvg');
end