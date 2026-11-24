%% =======================================================================
% Figure S6C and Figure S7A, D, E, G, H, I
%
% Parameters of the best-fitting model (model 9) of each participant
%% =======================================================================

clear all
close all


%% =======================================================================
% PATHS
%% =======================================================================

% Folder where this script is located
scriptFolder = fileparts(mfilename('fullpath'));
if isempty(scriptFolder), scriptFolder = pwd; end

% Go up one folder (repository root)
projectFolder = fileparts(scriptFolder);

% Path to model results
modelFolder = fullfile(projectFolder, 'data', 'model', 'results');


%% =======================================================================
% LOAD THE PARAMETERS AND COMPUTE ADAPTATION KERNEL AND DECISION BOUND
%% =======================================================================

nSubjects = 20;

allParameters  = nan(nSubjects,11);
adaptation     = nan(nSubjects,10);
decisionBound  = nan(nSubjects,10);

for subject = 1:nSubjects

    load(fullfile(modelFolder, ['Parameters_model_9_sbj_' num2str(subject) '.mat']));
    parameters = Parameters(subject,:,2);

    % Noise in the units of the model
    parameters(2) = parameters(2) / 20;

    allParameters(subject,:) = parameters;

    % Adaptation kernel (normalized to sum 1), as in model_simulation.m
    tau = parameters(3) * 15 + 0.0001;
    adaptation(subject,:) = parameters(4) + (1 - parameters(4)) .* exp(-(0:9) ./ tau);
    adaptation(subject,:) = adaptation(subject,:) ./ sum(adaptation(subject,:));

    % Collapsing decision bound, as in model_simulation.m
    initialBound  = abs(parameters(8)) * 15;
    collapseRate  = abs(parameters(6) * 50);
    decisionBound(subject,:) = initialBound - initialBound .* ((0:9) ./ ((0:9) + 1/collapseRate));

end


%% =======================================================================
% FIGURE 1 (Figure S7A) - ADAPTATION KERNEL
%% =======================================================================

figure(1)
hold on

plot(1:10, adaptation', 'color', [0.7 0.7 0.7], 'linewidth', 1);
plot(1:10, mean(adaptation,1), 'color', 'k', 'linewidth', 5);

title(sprintf('Average \\tau: %.2f', mean(allParameters(:,3) * 15)), 'FontWeight', 'normal')
xticks([1 10])
xlabel('Sample number')
ylabel('Adapted value (a.u.)')
box off
set(gca, 'fontsize', 15)


%% =======================================================================
% FIGURE 2 (Figure S7G) - DECISION BOUND
%% =======================================================================

figure(2)
hold on

plot(1:10, decisionBound', 'color', [0.7 0.7 0.7], 'linewidth', 1);
plot(1:10, mean(decisionBound,1), 'color', 'k', 'linewidth', 5);

xticks([1 10])
xlabel('Sample number')
ylabel('Decision bound (a.u.)')
box off
set(gca, 'fontsize', 15)


%% =======================================================================
% FIGURE 3 (Figure S6C, S7D, E, H, I) - PARAMETER ESTIMATES
%% =======================================================================

parameterNumbers = [5, 11, 1, 2, 10];
parameterNames   = {'Feed-forward inhibition (S6C)', 'Bias (S7D)', 'Leak (S7E)', ...
                    'Noise (S7H)', 'Confidence asymmetry weight (S7I)'};

figure(3)

for panel = 1:numel(parameterNumbers)

    subplot(1,numel(parameterNumbers),panel)
    hold on

    values = allParameters(:,parameterNumbers(panel));

    bar(1, mean(values), 'facecolor', [0.96 0.60 0.60])

    scatter(1 + (rand(nSubjects,1)*2 - 1)*0.05, values, 40, 'k', 'filled', ...
        'markerfacealpha', 0.5)

    xticks([])
    ylabel('Parameter estimate (a.u.)')
    title(parameterNames{panel}, 'FontWeight', 'normal', 'FontSize', 10)
    set(gca, 'fontsize', 12)

end
