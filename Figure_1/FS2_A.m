%% =======================================================================
% Figure S2A
%
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

% Path to functions (including the VBA toolbox)
addpath(genpath(fullfile(projectFolder, 'functions')));


%% =======================================================================
% SETTINGS
%% =======================================================================

nModels   = 19;
nSubjects = 20;

logLikelihood = nan(nSubjects,nModels);


%% =======================================================================
% LOG-LIKELIHOOD OF EACH MODEL AND SUBJECT
%% =======================================================================

for model = 1:nModels
    for subject = 1:nSubjects

        load(fullfile(modelFolder, ...
            ['All_data_model_' num2str(model) '_sbj_' num2str(subject) '.mat']));

        % Observed report of each trial (trials x 4, one 1 per row)
        % 1 = right/high, 2 = right/low, 3 = left/high, 4 = left/low
        observedReport = zeros(size(All_Conf));
        observedReport(Real_data(:,1)==1 & Real_data(:,2)==2,1) = 1;
        observedReport(Real_data(:,1)==1 & Real_data(:,2)==1,2) = 1;
        observedReport(Real_data(:,1)==0 & Real_data(:,2)==2,3) = 1;
        observedReport(Real_data(:,1)==0 & Real_data(:,2)==1,4) = 1;

        % Model probabilities of the 4 reports (no zero probabilities)
        modelProbability = double(All_Conf) * 0.99 + 0.005;

        % Remove trials without a report
        noReport = sum(observedReport,2) == 0;
        observedReport(noReport,:)   = [];
        modelProbability(noReport,:) = [];

        % Average log-likelihood per trial
        logLikelihood(subject,model) = ...
            nanmean(log(nansum(observedReport .* modelProbability, 2)));

    end
end


%% =======================================================================
% RANDOM-EFFECTS BAYESIAN MODEL COMPARISON
%% =======================================================================

options = [];
options.DisplayWin = 0;
options.verbose    = 0;

[~, out] = VBA_groupBMC(logLikelihood', options);


%% =======================================================================
% FIGURE (Figure S2A) - EXCEEDANCE PROBABILITY
%% =======================================================================

figure
hold on

bar(out.ep, ...
    'FaceColor', 0.75 * ones(1,3), ...
    'BarWidth', 0.6, ...
    'EdgeColor', 0.5 * ones(1,3), ...
    'LineWidth', 1.5);

ylim([0 1.1])
set(gca, ...
    'XTick', 1:nModels, ...
    'YTick', [0 0.25 0.5 0.75 1], ...
    'FontSize', 15, ...
    'Box', 'off', ...
    'TickDir', 'out');
ylabel('Exceedance probability')
xlabel('Model')
