%experiment6_eve_hmm: Simulate an HMM session and infer the parameters of
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
file_path = 'data/experiment6_eve_hmm';

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
          FontSize, FontName, Interpreter, FigureWidth, ...
          method, display, chunkSize};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (~loadData)
    % Simulate the session
    C = simulate(N, true_thetas{:});

    % Assume HMM
    fprintf('\nAssuming HMM...\n');
    [thetaP, varR, thetaR, thetaR_MAP, samples, Ks, K_theory] = post_process(C, params);

    % Assume iid
    fprintf('\nAssuming i.i.d. ...\n');

    params_iid = params;
    params_iid{1}{2}{1} = 0;
    params_iid{1}{2}{2} = 0;
    
    [thetaP_iid, varR_iid, thetaR_iid, thetaR_MAP_iid, samples_iid, Ks_iid] = post_process(C, params_iid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save/Load Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    close all; % close all open windows
    
    if (~isempty(file_path))
        
        fprintf('Saving data ...\n');

        save(file_path, ...
             'samples', 'samples_iid',...
             'varR', 'varR_iid',...
             'thetaR', 'thetaR_iid', ...
             'thetaP', 'thetaP_iid', ...
             'thetaR_MAP', 'thetaR_MAP_iid', ...
             'Ks', 'Ks_iid', ...
             'K_theory');
    end

else
    fprintf('Loading data ...\n');
    load(file_path);

    fprintf('\nAssuming HMM...\n');
    printResults(thetaR_MAP, thetaR, varR);

    fprintf('\nAssuming i.i.d. ...\n');
    printResults(thetaR_MAP_iid, thetaR_iid, varR_iid);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get LaTeX-formatted labels
[param_labels, K_labels] = getLabels(Nl, varR, Interpreter);
param_labels_iid = getLabels(Nl, varR_iid, Interpreter);

% Display posterior distributions of parameters assuming HMM
displayHistograms(samples, thetaP{3}, thetaP{4}, ...
                  'ground_truth', thetaR, ...
                  'labels', param_labels, ...
                  'CIp', CIp, ...
                  'FontSize', FontSize, ...
                  'FontName', FontName, ...
                  'Interpreter', Interpreter, ...
                  'FigureWidth', FigureWidth);

set(gcf, 'Name', 'Parameters - HMM', 'NumberTitle', 'off');
eve_hmm_params_plot = gcf;

% Display posterior distributions of parameters assuming i.i.d.
displayHistograms(samples_iid, thetaP_iid{3}, thetaP_iid{4}, ...
                  'ground_truth', thetaR_iid, ...
                  'labels', param_labels_iid, ...
                  'CIp', CIp, ...
                  'FontSize', FontSize, ...
                  'FontName', FontName, ...
                  'Interpreter', Interpreter, ...
                  'FigureWidth', FigureWidth);

set(gcf, 'Name', 'Parameters - i.i.d.', 'NumberTitle', 'off');
eve_iid_params_plot = gcf;

% Display distributions of secure key rates for each intensity assuming HMM
displayHistograms(Ks, ones(1, Nl), zeros(1, Nl), ...
                  'ground_truth', K_theory, ...
                  'labels', K_labels, ...
                  'CIp', CIp, ...
                  'FontSize', FontSize, ...
                  'FontName', FontName, ...
                  'Interpreter', Interpreter, ...
                  'FigureWidth', FigureWidth);

set(gcf, 'Name', 'Key Rates - HMM', 'NumberTitle', 'off');
eve_hmm_Ks_plot = gcf;

% Display distributions of secure key rates for each intensity assuming i.i.d.
displayHistograms(Ks_iid, ones(1, Nl), zeros(1, Nl), ...
                  'ground_truth', K_theory, ...
                  'labels', K_labels, ...
                  'CIp', CIp, ...
                  'FontSize', FontSize, ...
                  'FontName', FontName, ...
                  'Interpreter', Interpreter, ...
                  'FigureWidth', FigureWidth);

set(gcf, 'Name', 'Key Rates - i.i.d.', 'NumberTitle', 'off');
eve_iid_Ks_plot = gcf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save Figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if savePlots
    fprintf('Saving plots ...\n');
    print(eve_hmm_params_plot, 'figures/eve_hmm_params', '-dpdf');
    print(eve_iid_params_plot, 'figures/eve_hmm_params_with_iid', '-dpdf');

    print(eve_hmm_Ks_plot, 'figures/eve_hmm_Ks', '-dpdf');
    print(eve_iid_Ks_plot, 'figures/eve_hmm_Ks_with_iid', '-dpdf');
end
