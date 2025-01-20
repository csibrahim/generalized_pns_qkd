%experiment4_rates: Plots secure key rate vs. distance for the decoy, 
%                   corrected decoy, and proposed protocols across different
%                   intensity levels.
%
% Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
% Licensed under the MIT License (see LICENSE file for full details).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup                                                                  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Fix the random seed for reproducibility
seed = 1;

% Secure-key generation rate parameters
q = 0.5;                   % QKD efficiency
f = 1.22;                  % Error correction efficiency factor

savePlots = true;          % Set to 'true' to save the plot

rng(seed);                 % Set the random seed
restoredefaultpath;        % Restore default path and clear workspace
addpath(genpath('../.'));  % Add all subfolders to the search path

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

scale = 2;                    % Factor scaling for screen readability
FigureWidth = 180 ;           % Width of the figure in mm
FontSize = 8;                 % Font size for plot labels and text
FontName = 'Times New Roman'; % Font name for plot labels and text
aspectRatio = 2/5;            % Ratio of height/width
plot_spacing = 0.1;           % Spacing in km for theoretical distances
majorTick = 10;               % Major-ticks spacing
minorTick = 5;                % Minor-ticks spacing

% Scale the figure size and font
FigureWidth = FigureWidth * scale;
FontSize = FontSize * scale;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bob's Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Detection parameters for Bob's detectors (based on GYS protocol)
pas = [0 0.05 0.1];        % Array of after-pulsing probabilities
pc = 0.045;                % Detector efficiency
pd = 1.7e-6;               % Dark count probability
pe = 0.033;                % Misalignment probability

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Alice's Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha = 0.21;              % Attenuation coefficient
dAB = 220;                 % Maximum distance between Alice and Bob in km

% Decoy state intensities
mu = 0.48;                 % Signal intensity
nu1 = 0.05;                % First decoy state intensity
nu2 = 0;                   % Second decoy state intensity

% Proposed method additional intensities
mus = [1 5 10];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare Figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x_label = 'Distance in km';
y_label = 'K ';
prepFigure(FigureWidth, aspectRatio, FontSize, FontName, x_label, y_label);

lineWidth = FigureWidth / 200;

hold on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Theoretical Values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dABs = 0:plot_spacing:dAB; % Define distance range
Nl = length(mus);          % Number of additional intensities
colors = hsv(Nl);          % Distinct color for each lambda

% Decoy-state protocol plot
[dABs_decoy, Ks_decoy] = K_decoy(dABs, alpha, mu, nu1, nu2, pd, pe, pc, q, f);
h_decoy = plot(dABs_decoy, Ks_decoy, ':', ...
               'Color', 'k', ...
               'LineWidth', lineWidth);

% Corrected decoy-state protocol plot
[dABs_corrected, Ks_corrected] = K_corrected(dABs, alpha, mu, nu1, nu2, pd, pe, pc, q, f);
h_corrected = plot(dABs_corrected, Ks_corrected, '--', ...
                   'Color', 'k', ...
                   'LineWidth', lineWidth);

% Proposed protocol with mu
[dABs_proposed, Ks_proposed] = K_proposed(dABs, alpha, mu, pd, pe, pc, q, f);
h_mu = plot(dABs_proposed, Ks_proposed, '-', ...
            'Color', 'k', ...
            'LineWidth', lineWidth);

% Plots for additional intensities for the proposed protocol
h_lambdas = cell(1, Nl);
for i = 1:Nl
    [dABs_proposed, Ks_proposed] = K_proposed(dABs, alpha, mus(i), pd, pe, pc, q, f);
    h_lambdas{i} = plot(dABs_proposed, Ks_proposed, '-', ...
                        'Color', colors(i, :), ...
                        'LineWidth', lineWidth);
end

% Legend labels
intensities_legends = arrayfun(@(x) sprintf('Proposed (\\mu = %d)', x), mus, 'UniformOutput', false);
decoy_legend = {['Decoy (\mu = ', num2str(mu), ')']};
corrected_legend = {['Corrected Decoy (\mu = ', num2str(mu), ')']};
proposed_legend = {['Proposed (\mu = ', num2str(mu), ')']};

% Legend handles
legends = [decoy_legend, corrected_legend, proposed_legend, intensities_legends(:)'];
handles = [h_decoy, h_corrected, h_mu, h_lambdas{:}];

% Plot legends
legend(handles, legends, ...
       'Location', 'northoutside', ....
       'Orientation', 'horizontal', ....
       'FontSize', FontSize, ...
       'FontName', FontName, ...
       'NumColumns', 3);

% Set the y-axis scale to logarithmic
set(gca, 'YScale', 'log');

% Adjust x and y-axis limits
ax = gca;
lines = findall(ax, 'Type', 'Line');

minYs = arrayfun(@(h) min(h.YData), lines);
maxYs = arrayfun(@(h) max(h.YData), lines);

% Set the lower limit to 0
ax.YLim(1) = max(minYs);

% ylim([max(minYs) max(maxYs)]);
xlim([0 dAB]);

% Display and format axes grid
grid on;
ax.XMinorGrid = 'on';
ax.XMinorTick = 'on';

ax.XAxis.TickValues = 0:majorTick:dAB;
ax.YAxis.TickValues = 10.^(ceil(log10(max(minYs))):0);
ax.XAxis.MinorTickValues = minorTick:minorTick:dAB;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save Figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if savePlots
    fprintf('Saving the plot ...\n');
    print('figures/key_rates', '-dpdf');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Secure-key Rates Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dABs, Ks] = K_decoy(dABs, alpha, mu, nu1, nu2, pd, pe, pc, q, f)
    % K_decoy: Computes secure key rate for decoy-state protocol
    %          given distances and protocol parameters.

    e0 = 1 / 2; % Background error
    pABs = dist2prob(dABs, alpha);
    eta = pABs * pc;

    Qmu  = pd + 1 - exp(-eta * mu);
    Qnu1 = pd + 1 - exp(-eta * nu1);
    Qnu2 = pd + 1 - exp(-eta * nu2);

    EQmu = e0 * pd + pe * (1 - exp(-eta * mu));
    EQnu1 = e0 * pd + pe * (1 - exp(-eta * nu1));
    EQnu2 = e0 * pd + pe * (1 - exp(-eta * nu2));

    delta_mu = EQmu ./ Qmu;

    Y0 = max( (nu1 * Qnu2 * exp(nu2) - nu2 * Qnu1 * exp(nu1)) / (nu1 - nu2) ,0);
    
    Y1 = (mu / (mu * nu1 - mu * nu2 - nu1^2 + nu2^2)) * ...
         (Qnu1 * exp(nu1) - Qnu2 * exp(nu2) - (nu1^2 - nu2^2) * (Qmu*exp(mu) - Y0) / mu^2);
    
    Q1 = Y1 * mu * exp(-mu);
    e1 = (EQnu1 * exp(nu1) - EQnu2 * exp(nu2)) ./ ((nu1 - nu2) * Y1);

    Ks = q * (-f * Qmu .* H2(delta_mu) + Q1 .* (1 - H2(e1)));

    % Remove non-postive rates
    fail = Ks < 0 | isnan(Ks) | ~isreal(Ks);
    dABs(fail) = [];
    Ks(fail) = [];

end

function [dABs, Ks] = K_corrected(dABs, alpha, mu, nu1, nu2, pd, pe, pc, q, f)
    % K_corrected: Computes secure key rate for corrected decoy-state protocol
    %              given distances and protocol parameters.

    Ks = zeros(size(dABs));

    n = length(dABs);

    thetaB = {0, 0, pc, pc, pd, pd, pe};
    thetaE = {0, 0, 0, 0};

    for i=1:n

        thetaA = {mu, alpha, dABs(i)}; 
        [EQmu, Qmu] = EQ(thetaA, thetaB, thetaE);
        EQmu = 2 * EQmu(2);
        Qmu = 2 * Qmu(2);

        thetaA = {nu1, alpha, dABs(i)}; 
        [EQnu1, Qnu1] = EQ(thetaA, thetaB, thetaE);
        EQnu1 = 2 * EQnu1(2);
        Qnu1 = 2 * Qnu1(2);

        thetaA = {nu2, alpha, dABs(i)}; 
        [EQnu2, Qnu2] = EQ(thetaA, thetaB, thetaE);
        EQnu2 = 2 * EQnu2(2);
        Qnu2 = 2 * Qnu2(2);

        delta_mu = EQmu ./ Qmu;

        Y0_L = max( (nu1 * Qnu2 * exp(nu2) - nu2 * Qnu1 * exp(nu1)) / (nu1 - nu2) ,0);
        
        Y1 = (mu / (mu * nu1 - mu * nu2 - nu1^2 + nu2^2)) * ...
             (Qnu1 * exp(nu1) - Qnu2 * exp(nu2) - (nu1^2 - nu2^2) * (Qmu*exp(mu) - Y0_L) / mu^2);
        
        Q1 = Y1 * mu * exp(-mu);
        e1 = (EQnu1 * exp(nu1) - EQnu2 * exp(nu2)) ./ ((nu1 - nu2) * Y1);

        Ks(i) = q * (-f * Qmu .* H2(delta_mu) + Q1 .* (1 - H2(e1)));

    end

    % Remove non-postive rates
    fail = Ks < 0 | isnan(Ks) | ~isreal(Ks);
    dABs(fail) = [];
    Ks(fail) = [];

end

function [dABs, Ks] = K_proposed(dABs, alpha, mu, pd, pe, pc, q, f)
    % K_proposed: Computes secure key rate for the proposed protocol
    %             given distances and protocol parameters.
    
    Ks = zeros(size(dABs));

    n = length(dABs);

    for i=1:n

        % EQ = 1-(1-pd)*exp(-mu*pABs(i)*pc*pe);
        %  Q = 1-(1-pd)^2*exp(-mu*pABs(i)*pc);

        thetaA = {mu, alpha, dABs(i)}; 
        thetaB = {0, 0, pc, pc, pd, pd, pe};
        thetaE = {0, 0, 0, 0};
        
        [EQmu, Qmu] = EQ(thetaA, thetaB, thetaE, false);
        EQmu = 2 * EQmu(2);
        Qmu = 2 * Qmu(2);

        delta = EQmu ./ Qmu;
        Delta = 0;
        
        Ks(i) = q * Qmu .* (- f * H2(delta) + (1 - Delta) .* (1-H2(delta ./ (1 - Delta))));

    end
    
    % Remove non-postive rates
    fail = Ks < 0 | isnan(Ks) | ~isreal(Ks);
    dABs(fail) = [];
    Ks(fail) = [];
    
end