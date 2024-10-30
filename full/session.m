%session: This script sets up and simulates a quantum key distribution session, 
%         allowing estimation of the parameters associated with an eavesdropper 
%         using Bayesian inference and MCMC techniques. It saves and loads data, 
%         performs optimization, and computes secure key rates.
%
% Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
% Licensed under the MIT License (see LICENSE file for full details).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup                                                                  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng(1);                    % Set the random seed for reproducibility

restoredefaultpath;clear   % Restore default class path
addpath(genpath('.'));     % Add all subfolders to the search path

% File path for loading/saving results (can be empty [])
file_path = 'session'; 

% Set to 'true' to load data, 'false' to save after computation
loadData = false;

% Variables options
variables = {'lambdas','alpha','dAB',...
             'pa0','pa1','pc0','pc1','pd0','pd1','pe',...
             'dAE','pEB','k','Delta'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Session Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Pulse options
N = 1e8;          % Number of pulses in each simulation run
Nl = 9;           % Number of intensity levels (lambdas)

varR = {'dAE','pEB','k','Delta'};          % Random variables
varF = setdiff(variables, varR, 'stable'); % Fixed variables

% Get LaTeX-formatted labels
[param_labels, R_labels] = getLabels(Nl, varR);

% Variance options
noise = 0;  % noise to add during simulation
sigma = 0;  % standard deviation of the random variables

% Optimization options
maxEpochs = 50;   % Maximum epochs for outer optimization
maxIters  = 200;  % Maximum iterations for inner optimization

% Sampling options
Ns = 1e4;         % Number of samples
Nb = 1e3;         % Number of burn-in samples
display = false;  % Live histogram update (slow)
method = 'srss';  % Sampler options: 
                  %     'cmss'  - Covariance-Matching Slice Sampler
                  %     'srss'  - Shrinking-Rank Slice Sampler
                  %     'slice' - Slice sampling (within Gibbs)

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

alpha = 0.21; % Attenuation coefficient in channel
dAB = 100;    % Distance between Alice and Bob (in km)

% Get optimal intensity levels
lambdas = getLambdas(Nl, alpha, dAB, thetaB);

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
% Process Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Collect parameters for Alice, Bob, and Eve
thetas = {thetaA, thetaB, thetaE};

% Process random and fixed variables
[varR, varF] = processVars(thetas, varF, varR, sigma); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(~loadData)

    [C,D0,D1,l,a,b,x,e,sample_thetas] = simulate(N, thetas, varR, 'noise', noise);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process and get Prior Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [thetaP, thetaR, thetaF, varR, varF] = processParams(thetas, sample_thetas, ...
                                                            varF, varR, sigma);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Maximum a posteriori probability (MAP) estimation through optimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    thetaR_MAP = MAP(C, varR, varF, thetaF, thetaP, ...
                     'maxEpochs', maxEpochs, ...
                     'maxIters', maxIters);

    printResults(thetaR_MAP, thetaR, varR);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical Integration through MCMC sampling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    [samples,logPDFs] = sample(method, thetaR_MAP, Ns, Nb, ...
                        C, varR, varF, thetaF, thetaP, ...
                        'ground_truth', thetaR, ...
                        'display',display, ...
                        'labels', param_labels);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation of Secure Key Rate and Display Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Calculate the secure-key rate, gain and QBER given the true parameters
    [delta, Q, Delta] = deltaQ(thetaR, thetaF, varR, varF);
    R_theory = R(Delta, delta, Q);
    
    % Calculate the secure-key rate, gain and QBER given the parameters of the samples
    [deltas, Qs, Deltas] = deltaQ(samples, thetaF, varR, varF);
    Rs = R(Deltas, deltas, Qs);
    
    % Discard the non-matching basis
    R_theory = R_theory(Nl+1:end);
    Rs = Rs(:, Nl+1:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save/Load Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    close all;
    
    if(~isempty(file_path))
        save(file_path, 'samples', 'thetaR', 'thetaP', 'R_theory', 'Rs');
    end

else
    load(file_path);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Display the posterior distributions of the parameters
[upper_params, medians_params, lowers_params] = displayHistograms(samples, thetaP{3}, thetaP{4}, ...
                                                                  'ground_truth', thetaR, ...
                                                                  'labels', param_labels);

% Display the distributions of the secure key rates for each intensity
[upper_Rs, median_Rs, lower_Rs] = displayHistograms(Rs, ones(1,Nl), zeros(1,Nl), ...
                                                    'ground_truth', R_theory, ...
                                                    'labels', R_labels);