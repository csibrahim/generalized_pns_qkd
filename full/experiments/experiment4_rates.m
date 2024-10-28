% experiment4_rates: Plots secure key rate vs. distance for various protocols 
%                    including decoy, corrected decoy, and proposed protocols 
%                    across different intensity levels.
%
% Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
% Licensed under the MIT License (see LICENSE file for full details).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup                                                                  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng(1);                    % Set the random seed for reproducibility

restoredefaultpath; clear; % Restore default class path and clear workspace
addpath(genpath('../.'));  % Add all subfolders to the search path

% Secure-key generation rate parameters
q = 0.5;  % QKD efficiency
f = 1.22; % Error correction efficiency factor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FontSize = 30;                   % Font size for plot labels and text
FigureWidth = 1035;              % Width of the figure in points
plot_spacing = 0.1;              % Spacing in km for theoretical distance calculations 

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
dAB = 250;                       % Maximum distance between Alice and Bob in km

% Decoy state intensities
mu = 0.48;                       % Signal intensity
nu1 = 0.05;                      % First decoy state intensity
nu2 = 0;                         % Second decoy state intensity

% Proposed method additional intensities
lambdas = [1 5 10];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Theoretical Values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

prepFigure(FigureWidth, FontSize, 'Distance in km', 'R');

legendFontSize = FontSize * 0.75;
lineWidth = FigureWidth / 400;

% Define distance range
dABs = 0:plot_spacing:dAB;
hold on;

% Decoy-state protocol plot
[dABs_decoy, Rs_decoy] = R_decoy(dABs, alpha, mu, nu1, nu2, pd, pe, pc, q, f);
h_decoy = plot(dABs_decoy, Rs_decoy, 'k:', 'LineWidth', lineWidth);

% Corrected decoy-state protocol plot
[dABs_corrected, Rs_corrected] = R_corrected(dABs, alpha, mu, nu1, nu2, pd, pe, pc, q, f);
h_corrected = plot(dABs_corrected, Rs_corrected, 'k--', 'LineWidth', lineWidth);

% Proposed protocol with mu
[dABs_proposed, Rs_proposed] = R_proposed(dABs, alpha, mu, pd, pe, pc, q, f);
h_mu = plot(dABs_proposed, Rs_proposed, 'k-', 'LineWidth', lineWidth);

% Plot for each lambda intensity in proposed protocol
Nl = length(lambdas);
colors = hsv(Nl); % Distinct color for each lambda

h_lambdas = cell(1, Nl);
for i = 1:Nl
    [dABs_proposed, Rs_proposed] = R_proposed(dABs, alpha, lambdas(i), pd, pe, pc, q, f);
    h_lambdas{i} = plot(dABs_proposed, Rs_proposed, '-', ...
                        'LineWidth', lineWidth, ...
                        'Color', colors(i, :));
end

% Create legends for each protocol and intensity
intensities_legend = arrayfun(@(x) sprintf('Proposed ($\\mu=$%d)', x), lambdas, 'UniformOutput', false);
decoy_legend = {['Decoy ($\mu=', num2str(mu), '$)']};
corrected_legend = {['Corrected Decoy ($\mu=', num2str(mu), '$)']};
proposed_legend = {['Proposed ($\mu=', num2str(mu), '$)']};

legends = [decoy_legend, corrected_legend, proposed_legend, intensities_legend(:)'];
handles = [h_decoy, h_corrected, h_mu, h_lambdas{:}];
legend(handles, legends, ...
       'Location', 'northoutside', ....
       'Orientation', 'horizontal', ....
       'Interpreter', 'latex', ...
       'FontSize', legendFontSize, ...
       'NumColumns', 3);

set(gca, 'YScale', 'log');

% Adjust x and y-axis limits
fig = gcf;
axesHandle = fig.CurrentAxes;
lines = findall(axesHandle, 'Type', 'Line');

minYs = arrayfun(@(h) min(h.YData), lines);
maxYs = arrayfun(@(h) max(h.YData), lines);
maxXs = arrayfun(@(h) max(h.XData), lines);

max_y = 10^ceil(log10(max(maxYs)));
ylim([max(minYs) max_y]);
xlim([0 max(maxXs) * 1.05]);

function [dABs, Rs] = R_decoy(dABs, alpha, mu, nu1, nu2, pd, pe, pc, q, f)
    % R_decoy: Computes secure key rate for decoy-state protocol
    %           given distances and protocol parameters.

    e0 = 1 / 2; % Background error
    pABs = dist2prob(dABs, alpha);
    eta = pABs * pc;

    Qmu  = pd + 1 - exp(-eta * mu);
    Qnu1 = pd + 1 - exp(-eta * nu1);
    Qnu2 = pd + 1 - exp(-eta * nu2);

    EmuQmu = e0 * pd + pe * (1 - exp(-eta * mu));
    Enu1Qnu1 = e0 * pd + pe * (1 - exp(-eta * nu1));
    Enu2Qnu2 = e0 * pd + pe * (1 - exp(-eta * nu2));

    Emu = EmuQmu ./ Qmu;
    
    Y1 = (mu / (mu * nu1 - mu * nu2 - nu1^2 + nu2^2)) * ...
         (Qnu1 * exp(nu1) - Qnu2 * exp(nu2) - (nu1^2 - nu2^2) * (Qmu - pd) / mu^2);
    Q1 = Y1 * mu * exp(-mu);
    e1 = (Enu1Qnu1 * exp(nu1) - Enu2Qnu2 * exp(nu2)) ./ ((nu1 - nu2) * Y1);
    
    Rs = q * (-f * Qmu .* H2(Emu) + Q1 .* (1 - H2(e1)));

    % Remove negative rates
    dABs(Rs < 0) = [];
    Rs(Rs < 0) = [];

end

function [dABs, Rs] = R_corrected(dABs, alpha, mu, nu1, nu2, pd, pe, pc, q, f)
    % R_corrected: Computes secure key rate for corrected decoy-state protocol.

    pABs = dist2prob(dABs, alpha);
    eta = pABs * pc;

    Qmu  = 1 - power(1 - pd, 2) * exp(-eta * mu);
    Qnu1 = 1 - power(1 - pd, 2) * exp(-eta * nu1);
    Qnu2 = 1 - power(1 - pd, 2) * exp(-eta * nu2);

    EmuQmu = 1 - (1 - pd) * exp(-eta * mu * pe);
    Enu1Qnu1 = 1 - (1 - pd) * exp(-eta * nu1 * pe);
    Enu2Qnu2 = 1 - (1 - pd) * exp(-eta * nu2 * pe);

    Emu = EmuQmu ./ Qmu;
    
    Y1 = (mu / (mu * nu1 - mu * nu2 - nu1^2 + nu2^2)) * ...
         (Qnu1 * exp(nu1) - Qnu2 * exp(nu2) - (nu1^2 - nu2^2) * (Qmu - pd) / mu^2);
    Q1 = Y1 * mu * exp(-mu);
    e1 = (Enu1Qnu1 * exp(nu1) - Enu2Qnu2 * exp(nu2)) ./ ((nu1 - nu2) * Y1);
    
    Rs = q * (-f * Qmu .* H2(Emu) + Q1 .* (1 - H2(e1)));

    % Remove negative rates
    dABs(Rs < 0) = [];
    Rs(Rs < 0) = [];

end

function [dABs, Rs] = R_proposed(dABs, alpha, mu, pd, pe, pc, q, f)
    % R_proposed: Computes secure key rate for proposed protocol.
    
    varF = {'lambdas', 'alpha', ...
            'pa0', 'pa1', 'pc0', 'pc1', 'pd0', 'pd1', 'pe', ...
            'dAE', 'k', 'Delta'};
    varR = {'dAB', 'pEB'};

    pABs = dist2prob(dABs, alpha);

    thetaA = {mu, alpha}; 
    thetaE = {0, 0, 0};
    thetaR = {dABs', pABs'};
    thetaB = {0, 0, pc, pc, pd, pd, pe};
    thetaF = [thetaA, thetaB, thetaE];
    
    [deltaQ, Q] = deltaQ_iid(thetaR, thetaF, varR, varF, true);
    deltaQ = 2 * deltaQ(:, 2)';
    Q = 2 * Q(:, 2)';
    delta = deltaQ ./ Q;
    Delta = 0;
    
    Rs = q * Q .* (1 - Delta - f * H2(delta) - (1 - Delta) .* H2(delta ./ (1 - Delta)));
    
    % Remove negative rates
    dABs(Rs < 0) = [];
    Rs(Rs < 0) = [];
end