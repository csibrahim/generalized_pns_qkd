%experiment3_deltas: Plots the error probability across various distances 
%                    and intensity levels, highlighting both the expected 
%                    and observed detection error probabilities for a 
%                    range of after-pulsing probabilities compared to the
%                    decoy protocol.
%
% Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
% Licensed under the MIT License (see LICENSE file for full details).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup                                                                  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Fix the random seed for reproducibility
seed = 1;

% File path for loading/saving results (can be empty [])
file_path = 'data/experiment3_deltas'; 

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
aspectRatio = 1/3;   % Ratio of height/width
plot_spacing = 1;    % Spacing in km for theoretical distance calculations   
sample_spacing = 5;  % Spacing in km for simulation data points

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Session Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 1e9; % Number of pulses in each simulation run

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bob's Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Detection parameters for Bob's detectors (based on GYS protocol)
pas = [0 0.05 0.1];  % Array of after-pulsing probabilities
pc = 0.045;          % Detector efficiency
pd = 1.7e-6;         % Dark count probability
pe = 0.033;          % Misalignment probability

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Alice's Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mu = 0.48;          % Signal intensity
alpha = 0.21;       % Attenuation coefficient
dAB = 150;          % Maximum distance between Alice and Bob in km

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare Figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x_label = 'Distance in km';
y_label = '$\delta$ ';

prepFigure(FigureWidth, aspectRatio, FontSize, x_label, y_label);

lineWidth = FigureWidth / 400;
markerSize = FigureWidth / 150;

hold on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Theoretical Values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Distances and transmission probabilities for theoretical calculations
dABs = (0:plot_spacing:dAB)';
pABs = dist2prob(dABs, alpha);

% Number of after-pulse scenarios
n_pa = length(pas); 

% Initialize cells to store theoretical error rates for each after-pulsing 
% probability and the decoy-state
Qs = cell(1, n_pa + 1);
deltas = cell(1, n_pa + 1);

% Theoretical error rates for the decoy protocol without after-pulsing
e0 = 1 / 2; % Background error probability 
Qs{1} = 0.5 * (pd + 1 - exp(-mu * pABs * pc));
deltas{1} = 0.5 * (e0 * pd + pe * (1 - exp(-mu * pABs * pc ))) ./ Qs{1};

% Correct theoretical error rates without after-pulsing
Qs{2} = 0.5 * (1-(1-pd)^2 * exp(-mu * pABs * pc));
deltas{2} = 0.5 * (1-(1-pd) * exp(-mu * pABs * pc * pe)) ./ Qs{2};

% Calculate theoretical probabilities for each after-pulsing probability setting
for i = 2:n_pa

    Qs{i + 1} = zeros(size(Qs{1}));
    deltas{i + 1} = zeros(size(deltas{1}));

    for j=1:length(dABs)

        thetaA = {mu, alpha, dABs(j)};
        thetaB = {pas(i), pas(i), pc, pc, pd, pd, pe};
        thetaE = {dABs(j), 1, 0, 0};
        [EQij, Qij] = EQ(thetaA, thetaB, thetaE, true);

        Qs{i + 1}(j) = Qij(2);
        deltas{i + 1}(j) = EQij(2) ./ Qij(2);

    end
    
end

% Plot confidence intervals and expected values for the decoy protocol and each after-pulsing probability
for i = 1:numel(Qs)
    
    E_S = N*Qs{i}; % The expected value of the signal probability

    % Expected value and variance of the error rate 
    E_rho = deltas{i};                           
    V_rho = (deltas{i}.*(1-deltas{i}))./E_S;

    % Beta distribution parameters with mean = E_rho and variance = V_rho
    [alphas, betas] = beta_parameters(E_rho, V_rho);

    % Copmute bonds from the parameters
    epsilon = (1-CIp) / 2;

    lower = betainv(epsilon,alphas,betas);    % Upper bound of confidence interval
    upper = betainv(1-epsilon,alphas,betas);  % Lower bound of confidence interval

    % Plot confidence interval as shaded regions
    h_CI = fill([dABs', fliplr(dABs')], [upper', fliplr(lower')], [0 0 0], ...
                'EdgeColor', 'k', ...
                'EdgeAlpha', 0.5, ...
                'FaceAlpha', 0.15);

    % Plot expected value
    h_E = plot(dABs, E_rho, 'k-', ...
               'LineWidth', lineWidth);

    % Make the line half transparent
    h_E.Color = [0 0 0 0.5]; 

    % Place a label at the end of each line
    if i == 1
        lineLabel = {'Decoy', '($p_{a}$=\,0\%)'};
    else
        lineLabel = ['$p_{a}$=\,', num2str(100 * pas(i - 1)), '\%'];
    end

    text(max(dABs) * 1.01, E_rho(end), lineLabel, ...
         'Interpreter', 'latex', ...
         'FontSize', FontSize);
end

% Define colors for each after-pulsing probability and initialize plot handles
colors = hsv(n_pa);
h_data = cell(1, n_pa);

for j = 1:n_pa
    h_data{j} = plot(nan, nan, 'o', ...
                     'Color', 'k', ...
                     'MarkerFaceColor', colors(j, :), ...
                     'MarkerSize', markerSize);
end

% Configure legend labels
CI_legend = {['$\mathrm{CI}_{', num2str(100*CIp), '\%}$']};
E_legend = {'Expected Value'};
data_legends = arrayfun(@(x) sprintf('$%d\\%%$', x), 100 * pas, 'UniformOutput', false);

% Configure legends handles
handles = [h_CI, h_E, h_data{:}];
legends = [CI_legend, E_legend, data_legends(:)'];

% Plot legends
legend(handles, legends, ...
       'Location', 'northoutside', ...
       'Orientation', 'horizontal', ...
       'Interpreter', 'latex', ...
       'FontSize', FontSize);

% Remove padding around the axes, but leave 10% to the east for legends
set(gca, 'LooseInset', [0, 0, 0.1, 0]);

% Get current axes handle
ax = gca;

% Set the lower limit to 0
ax.YLim(1) = 0;

% Convert y-axis to percentage and add '%'
ax.YAxis.TickLabels = strcat(string(yticks * 100), '\%');

% Setup minor ticks
ax.XAxis.MinorTickValues = 0:0.5*sample_spacing:dAB;
ax.YAxis.MinorTickValues = 0.01:0.02:max(yticks);

% A major x-tick every other sample
ax.XAxis.TickValues = 0:sample_spacing:dAB;

% Display and format axes grid
grid on;
ax.XMinorTick = 'on';
ax.XMinorGrid = 'on';
ax.YMinorTick = 'on';
ax.YMinorGrid = 'on';
ax.MinorGridLineStyle = ':';

% Force displaying the figure now
drawnow;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Measured Data Points (Simulation Results)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (loadData)
    % Load previously simulated data if specified
    fprintf('Loading data ...\n');
    load(file_path, 'dABs', 'measured');

    % Assign loaded data to handles
    for j = 1:n_pa
        set(h_data{j}, 'XData', dABs, 'YData', measured(:, j));
    end
else

    % Run simulation to generate new data if loading is not specified

    dABs = 0:sample_spacing:dAB;        % Disatnces to simulate at
    n_dABs = length(dABs);              % Number of distances
    measured = zeros(n_dABs, n_pa);     % Initialize array to store measurements

    Nl = 1; % Only one intensity (mu) is used for the key generation

    % Loop through each distance and after-pulsing setting
    for i = 1:length(dABs)
        for j = 1:n_pa
            
            % Print progress so far
            fprintf(['[', num2str((i - 1) * n_pa + j), '/', num2str(n_dABs * n_pa), ']']);
            
            % Update parameters for simulation
            thetaA = {mu, alpha, dABs(i)};
            thetaB = {pas(j), pas(j), pc, pc, pd, pd, pe};
            thetaE = {dABs(i), 1, 0, 0};

            % Run simulation
            [~, D0, D1, l, a, b, x] = simulate(N, thetaA, thetaB, thetaE);

            % Measure observed counts
            [R, S, B] = EQ_measure(Nl, D0, D1, l, a, b, x);

            % Consider double clicks for comparison with the decoy protocol
            R = R+B;
            S = S+B;

            rho = R ./ S;
            
            measured(i, j) = rho(2); % Consider only matching basis events

            % Plot the simulated data point
            set(h_data{j}, 'XData', dABs(1:i), 'YData', measured(1:i, j));

            % Update the figure
            drawnow;
            
            
        end
    end

    % Save simulated data
    if(~isempty(file_path))
        fprintf('Saving data ...\n');
        save(file_path, 'dABs', 'measured');
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save Figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if savePlots
    fprintf('Saving the plot ...\n');
    print('figures/error_rates', '-dpdf');
end