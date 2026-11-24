%% =======================================================================
% Figure S2B
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

% Fits of the synthetic datasets: rec_g09_d<dataset>_f<fitted model>.mat
recoveryFolder = fullfile(projectFolder, 'data', 'model', 'recovery');

%% =======================================================================
% SETTINGS
%% =======================================================================

generatingModel = 9;
fittedModels    = [3 6 9 11 14 18 19];

nModels = numel(fittedModels);


%% =======================================================================
% LOAD THE BIC OF EVERY DATASET AND FITTED MODEL
%% =======================================================================

files = dir(fullfile(recoveryFolder, ['rec_g' sprintf('%02d',generatingModel) '_d*_f*.mat']));

datasetNumber = nan(numel(files),1);
fittedNumber  = nan(numel(files),1);
BIC_file      = nan(numel(files),1);

for k = 1:numel(files)

    numbers = sscanf(files(k).name, ['rec_g' sprintf('%02d',generatingModel) '_d%d_f%d.mat']);

    datasetNumber(k) = numbers(1);
    fittedNumber(k)  = numbers(2);

    fit = load(fullfile(recoveryFolder, files(k).name), 'BIC');
    BIC_file(k) = fit.BIC;

end

% Keep the datasets fitted with all 7 models
datasets = unique(datasetNumber)';
complete = false(size(datasets));

for k = 1:numel(datasets)
    complete(k) = all(ismember(fittedModels, fittedNumber(datasetNumber == datasets(k))));
end

datasets  = datasets(complete);
nDatasets = numel(datasets);

% BIC: datasets x fitted models
BIC = nan(nDatasets,nModels);

for k = 1:nDatasets
    for model = 1:nModels
        BIC(k,model) = BIC_file(datasetNumber == datasets(k) & fittedNumber == fittedModels(model));
    end
end


%% =======================================================================
% WINNING MODEL AND BIC DIFFERENCE TO MODEL 9
%% =======================================================================

[~, winner] = min(BIC, [], 2);
nWins = accumarray(winner, 1, [nModels 1])';

deltaBIC      = BIC - BIC(:, fittedModels == generatingModel);
deltaBIC_mean = mean(deltaBIC,1);
deltaBIC_sem  = std(deltaBIC,0,1) ./ sqrt(nDatasets);


%% =======================================================================
% FIGURE (Figure S2B)
%% =======================================================================

modelNames = strcat('Model', {' '}, arrayfun(@num2str, fittedModels, 'UniformOutput', false));

barColors = repmat([0.75 0.75 0.75], nModels, 1);
barColors(fittedModels == generatingModel,:) = [0.85 0.2 0.2];

figure('Color', 'w', 'Position', [100 100 760 380]);

% ---------------------------------------------------------------
% Left: number of datasets won
% ---------------------------------------------------------------

subplot(1,2,1)
hold on

for model = 1:nModels

    barh(model, nWins(model), 0.6, 'FaceColor', barColors(model,:), 'EdgeColor', 'none');

    text(nWins(model) + 1, model, sprintf('%d/%d', nWins(model), nDatasets), ...
        'VerticalAlignment', 'middle', 'FontSize', 10);

end

set(gca, 'YTick', 1:nModels, 'YTickLabel', modelNames, 'YDir', 'reverse', ...
    'XAxisLocation', 'top', 'Box', 'off', 'TickDir', 'out', 'FontSize', 11);
ylim([0.4 nModels+0.6])
xlim([0 100])
xlabel('Number of datasets won')

% ---------------------------------------------------------------
% Right: BIC difference to model 9
% ---------------------------------------------------------------

subplot(1,2,2)
hold on

barh(1:nModels, deltaBIC_mean, 0.6, 'FaceColor', [0.6 0.6 0.6], 'EdgeColor', 'none');

errorbar(deltaBIC_mean, 1:nModels, deltaBIC_sem, 'horizontal', ...
    'LineStyle', 'none', 'Color', 'k', 'LineWidth', 1);

set(gca, 'YTick', 1:nModels, 'YTickLabel', [], 'YDir', 'reverse', ...
    'XAxisLocation', 'top', 'Box', 'off', 'TickDir', 'out', 'FontSize', 11);
ylim([0.4 nModels+0.6])
xlabel('\DeltaBIC (model - model 9)')
