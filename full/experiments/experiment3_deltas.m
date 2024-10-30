% experiment3_deltas: Plots the error probability across various distances and 
%                     intensity levels, highlighting both the expected and observed 
%                     detection error probabilities for a range of after-pulsing
%                     probabilities compared to the decoy protocol.
%
% Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
% Licensed under the MIT License (see LICENSE file for full details).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup                                                                  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng(1);                    % Set the random seed for reproducibility

restoredefaultpath; clear; % Restore default class path and clear workspace
addpath(genpath('../.'));  % Add all subfolders to the search path

% File path for loading/saving results
file_path = 'saved_data/experiment3_deltas';

% Set to 'true' to load data, 'false' to save after computation
loadData = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Session Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 1e9;  % Number of pulses in each session simulation run

% Define variable sets for the theoretical model
varF = {'lambdas', 'alpha', ...
        'pa0', 'pa1', 'pc0', 'pc1', 'pd0', 'pd1', 'pe', ...
        'dAE', 'k', 'Delta'}; % Fixed variables
varR = {'dAB', 'pEB'};        % Random variables

% Variance options
noise = 0;  % noise to add during simulation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FontSize = 30;       % Font size for plot labels and text
FigureWidth = 1035;  % Width of the figure in points
width = 3;           % Width of the confidence interval in standard deviations
plot_spacing = 1;    % Spacing in km for theoretical distance calculations   
sample_spacing = 5;  % Spacing in km for simulation data points

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bob's Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Detection parameters for Bob's detectors (based on GYS protocol)
pas = [0 0.05 0.1];              % Array of after-pulsing probabilities
pc = 0.045;                      % Detector efficiency
pd = 1.7e-6;                     % Dark count probability
pe = 0.033;                      % Misalignment probability

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Alice's Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha = 0.21;                    % Attenuation coefficient in channel
dAB = 150;                       % Maximum distance between Alice and Bob in km

% Decoy state intensities for quantum key distribution
mu = 0.48;                       % Signal intensity
nu1 = 0.05;                      % First decoy state intensity
nu2 = 0;                         % Second decoy state intensity

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Theoretical Values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Prepare figure for plotting
prepFigure(FigureWidth, FontSize, 'Distance in km', '$\delta$');

legendFontSize = FontSize * 0.75;
lineWidth = FigureWidth / 400;
markerSize = FigureWidth / 100;

% Define distance range for theoretical calculations
dABs = (0:plot_spacing:dAB)';
pABs = dist2prob(dABs, alpha);   % Calculate probability of transmission over distance

% Define fixed parameter sets for Alice and Eve
thetaA = {mu, alpha};
thetaE = {0, 0, 0};
thetaR = {dABs, pABs};

hold on;
n_pa = length(pas);

% Initialize cells to store theoretical probabilities and errors for each after-pulsing probability
Qs = cell(1, n_pa + 1);
deltas = cell(1, n_pa + 1);

% Calculate theoretical error rates for the decoy protocol without after-pulsing
e0 = 1 / 2; % Background error probability 
Qs{1} = pd + 1 - exp(-pABs * pc * mu);
deltas{1} = (e0 * pd + pe * (1 - exp(-pABs * pc * mu))) ./ Qs{1};

% Calculate theoretical probabilities for each after-pulsing probability setting
for i = 1:n_pa
    thetaB = {pas(i), pas(i), pc, pc, pd, pd, pe};
    thetaF = [thetaA, thetaB, thetaE];

    [deltaQi, Qi] = deltaQ(thetaR, thetaF, varR, varF, true);
    
    Qs{i + 1} = 2 * Qi(:, 2);
    deltas{i + 1} = deltaQi(:, 2) ./ Qi(:, 2);
end

% Plot confidence intervals and expected values for the decoy protocol and each after-pulsing probability
for i = 1:numel(Qs)
    
    sigma = sqrt(V_error(N, Qs{i}, deltas{i})); % Calculate standard deviation for confidence interval
    upper = deltas{i} + width * sigma;          % Upper bound of confidence interval
    expected = deltas{i};                       % Expected error rate
    lower = deltas{i} - width * sigma;          % Lower bound of confidence interval

    % Plot confidence interval as shaded region
    h_CI = fill([dABs', fliplr(dABs')], [upper', fliplr(lower')], [0 0 0], ...
                'EdgeColor', 'k', 'EdgeAlpha', 0.15, 'FaceAlpha', 0.15);

    % Plot expected value
    h_E = plot(dABs, expected, 'k-', 'LineWidth', lineWidth);
    h_E.Color = [0 0 0 0.5];

    % Place a label at the end of each line
    if i == 1
        lineLabel = {'Decoy', '($p_{a}$=\,0\%)'};
    else
        lineLabel = ['$p_{a}$=\,', num2str(100 * pas(i - 1)), '\%'];
    end
    text(max(dABs) * 1.01, expected(end), lineLabel, 'Interpreter', 'latex', 'FontSize', legendFontSize);
end

drawnow;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Measured Data Points (Simulation Results)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define colors for each after-pulsing probability and initialize plot handles
colors = hsv(n_pa);
h_data = cell(1, n_pa);

if (loadData)
    % Load previously simulated data if specified
    load(file_path, 'dABs', 'measured');
    for j = 1:n_pa
        h_data{j} = plot(dABs, measured(:, j), 'o', 'Color', 'k', ...
                         'MarkerFaceColor', colors(j, :), 'MarkerSize', markerSize);
    end
else
    % Run simulation to generate new data if loading is not specified
    dABs = 0:sample_spacing:dAB;
    n_dABs = length(dABs);
    measured = zeros(n_dABs, n_pa);

    Nl = 1; % Only one intensity (mu) is used for the key generation
    
    % Loop through each distance and after-pulsing setting
    for i = 1:length(dABs)
        for j = 1:n_pa
            fprintf(['[', num2str((i - 1) * n_pa + j), '/', num2str(n_dABs * n_pa), ']']);
            
            % Update parameters for simulation
            thetaA = {mu, alpha, dABs(i)};
            thetaB = {pas(j), pas(j), pc, pc, pd, pd, pe};
            thetaE = {dABs(i), 1, 0, 0};
            thetas = {thetaA, thetaB, thetaE};

            % Run simulation and calculate observed probabilities
            [~, D0, D1, l, a, b, x] = simulate(N, thetas, varR, ...
                                              'match', true, ...
                                              'noise', noise);

            [m, M, both] = measure(Nl, D0, D1, l, a, b, x);

            % Consider double clicks for comparison with the decoy protocol
            m = m+both;
            M = M+both;

            delta = m ./ M;
            
            measured(i, j) = delta(2); % Consider only matching basis events

            % Plot simulated data points
            h_data{j} = plot(dABs(i), measured(i, j), 'o', 'Color', 'k', ...
                             'MarkerFaceColor', colors(j, :), 'MarkerSize', markerSize);

            drawnow;
        end
    end
    save(file_path, 'dABs', 'measured'); % Save simulated data
end

% Configure legend for the plot
data_legends = arrayfun(@(x) sprintf('$%d\\%%$', x), 100 * pas, 'UniformOutput', false);
CI_legend = {['$\pm$\,', num2str(width), '$\,\sigma$']};
E_legend = {'Expected Value'};

legends = [CI_legend, E_legend, data_legends(:)'];
handles = [h_CI, h_E, h_data{:}];

legend(handles, legends, 'Location', 'northoutside', 'Orientation', 'horizontal', ...
       'Interpreter', 'latex', 'FontSize', legendFontSize);