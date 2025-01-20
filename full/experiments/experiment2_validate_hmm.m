%experiment2_validate_hmm: Simulates detection events across various intensities 
%                          and compares observed distributions of clicks and 
%                          errors with theoretical distributions derived from
%                          both the i.i.d. and HMM model.
%
% Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
% Licensed under the MIT License (see LICENSE file for full details).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup                                                                  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Fix the random seed for reproducibility
seed = 1;

% File path for loading/saving results (can be empty [])
file_path = 'data/experiment2_validate_hmm'; 

loadData = false;    % Set to 'true' to load data, 'false' to save after computation
savePlots = true;    % Set to 'true' to save plots

rng(seed);                 % Set the random seed
restoredefaultpath;        % Restore default path and clear workspace
addpath(genpath('../.'));  % Add all subfolders to the search path

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting Options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

scale = 1.5;                    % Factor scaling for screen readability
FigureWidth = 180;              % Width of the figure in mm
FontSize = 8;                   % Font size for plot labels and text
FontName = 'Times New Roman';   % Font name for plot labels and text
CIp = 0.99;                     % Confidence interval threshold

% Scale the figure size and font
FigureWidth = FigureWidth * scale;
FontSize = FontSize * scale;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Session Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Variable settings
variables = {'lambdas','alpha','dAB', ...
             'pa0','pa1','pc0','pc1','pd0','pd1','pe',...
             'dAE','pEB','k','Delta'};

varR = {'dAE','pEB','k','Delta'};          % Random variables
varF = setdiff(variables, varR, 'stable'); % Fixed variables

% Pulse options
N  = 1e4;   % Number of pulses in each simulation run
Nr = 1e4;   % Total number of simulation runs
Nl = 2;     % Number of intensity levels (lambdas)
limit = 10; % Maximum intensity

% Variance options
noise = 0;  % Noise level for simulation
sigma = 0;  % Standard deviation of the random variables

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
% Process Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Collect parameters for Alice, Bob, and Eve
thetas = {thetaA, thetaB, thetaE};

% Process random and fixed variables
[varR, varF] = processVars(thetas, varF, varR); 

% Process parameter priors and obtain random variables
[thetaP, thetaR, thetaF, varR, varF] = processParams(thetas, thetas, varF, varR, sigma);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(~loadData)

    % Initialize  matrices for counts, errors, signals, number of intensities
    C = zeros(Nr, 8 * Nl);    % Clicks per intensity/basis-match/detection-event 
    R = zeros(Nr, 2 * Nl);    % Erroneous clicks per intensity/basis-match
    S = zeros(Nr, 2 * Nl);    % XOR clicks per intensity/basis-match
    B = zeros(Nr, 2 * Nl);    % AND clicks per intensity/basis-match
    
    % Run simulation across the specified number of runs
    for i = 1:Nr
        
        % Simulate pulses and measure detection events
        [C(i,:), D0, D1, l, a, b, x] = simulate(N, thetas{:}, 'print', false);

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

% Compute HMM detection probabilities and plot results
Ps_hmm = P_hmm(thetaR, thetaF, varR, varF);

% Compute i.i.d. detection probabilities and plot results
Ps_iid = P_iid(thetaR, thetaF, varR, varF);

plotPs(Ps_hmm, C, ...
         'FigureWidth', FigureWidth, ...
         'FontSize', FontSize, ...
         'FontName', FontName, ...
         'CIp', CIp, ...
         'Ps_ref', Ps_iid, ...
         'LegendName', {'HMM', 'i.i.d'});

set(gcf, 'Name', 'Probabilities - HMM vs i.i.d.', 'NumberTitle', 'off');
Ps_plot = gcf;


% Compute HMM error rates and gain, and visualize them
[EQs_hmm, Qs_hmm] = EQ_hmm(thetaR, thetaF, varR, varF);

% Compute i.i.d. error rates and gain, and visualize them
[EQs_iid, Qs_iid] = EQ_iid(thetaR, thetaF, varR, varF);

plotEQ(N, R, S, EQs_hmm, Qs_hmm, ...
         'FigureWidth', FigureWidth, ...
         'FontSize', FontSize, ...
         'FontName', FontName, ...
         'CIp', CIp, ...
         'EQs_ref', EQs_iid, ...
         'Qs_ref', Qs_iid, ...
         'LegendName', {'HMM', 'i.i.d'});

set(gcf, 'Name', 'Gain/Error - HMM vs i.i.d.', 'NumberTitle', 'off');
EQs_plot = gcf;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save Figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if savePlots
    fprintf('Saving plots ...\n');
    print(Ps_plot, 'figures/validate_hmm_Ps', '-dpdf');
    print(EQs_plot, 'figures/validate_hmm_EQs', '-dpdf');
end