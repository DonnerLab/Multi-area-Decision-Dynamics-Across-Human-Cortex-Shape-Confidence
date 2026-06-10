%% =======================================================================
% Figure S2A: Model comparison using group BMC
%% =======================================================================

clearvars
close all

%% Paths

% Folder where this script is located
scriptFolder = fileparts(mfilename('fullpath'));

% Go up two folders
repoRoot = (fileparts(scriptFolder));

% Standard project paths
dataFolder      = fullfile(repoRoot, 'data');
modelFolder     = fullfile(dataFolder, 'model', 'results');
functionsFolder = fullfile(repoRoot, 'functions');

% Add functions (including VBA toolbox if stored here)
addpath(genpath(functionsFolder));

%% Settings

nModels   = 19;
nSubjects = 20;

LL_all = nan(nSubjects, nModels);

%% Compute log-likelihoods

for model = 1:nModels

    for s = 1:nSubjects

        load(fullfile(modelFolder, ...
            sprintf('All_data_model_%d_sbj_%d.mat', model, s)));

        % Create a trial × 4 matrix with each column corresponding to:
        % 1 = Right / High confidence
        % 2 = Right / Low confidence
        % 3 = Left  / High confidence
        % 4 = Left  / Low confidence

        Confidence_choices = zeros(size(All_Conf));

        Confidence_choices(Real_data(:,1)==1 & Real_data(:,2)==2,1) = 1;
        Confidence_choices(Real_data(:,1)==1 & Real_data(:,2)==1,2) = 1;
        Confidence_choices(Real_data(:,1)==0 & Real_data(:,2)==2,3) = 1;
        Confidence_choices(Real_data(:,1)==0 & Real_data(:,2)==1,4) = 1;

        % Probabilities from model fits
        Choices_P = double(All_Conf) * 0.99 + 0.005;

        % Remove invalid trials
        invalidTrials = sum(Confidence_choices,2) == 0;

        Confidence_choices(invalidTrials,:) = [];
        Choices_P(invalidTrials,:) = [];

        % Average log-likelihood
        LL_all(s,model) = nanmean( ...
            log(nansum(Confidence_choices .* Choices_P, 2)));
    end

end

%% Bayesian Model Comparison

opt = [];
opt.DisplayWin = 0;

[~, out] = VBA_groupBMC(LL_all', opt);

%% Plot

figure

bar(out.ep, ...
    'FaceColor', 0.75 * ones(1,3), ...
    'BarWidth', 0.6, ...
    'EdgeColor', 0.5 * ones(1,3), ...
    'LineWidth', 1.5);

hold on

ylim([0 1.1])

set(gca, ...
    'XTick', 1:nModels, ...
    'XTickLabel', 1:nModels, ...
    'XTickLabelRotation', 30, ...
    'FontSize', 6, ...
    'Box', 'off', ...
    'TickDir', 'out', ...
    'LineWidth', 1, ...
    'YTick', [0 1]);

ylabel('Expected probability')
xlabel('Model')

% Display model comparison metrics
disp('Protected exceedance probabilities:')
disp(out.pxp)

pxp = (1 - out.bor) * out.ep + out.bor / length(out.ep);

disp('Adjusted exceedance probabilities:')
disp(pxp)