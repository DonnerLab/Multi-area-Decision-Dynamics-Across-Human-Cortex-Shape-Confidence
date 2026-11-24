%% =======================================================================
% Figure S3E-F
%
% =======================================================================

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

% Path to model results (behavior + model predictions of each subject)
modelFolder = fullfile(projectFolder, 'data', 'model', 'results');

% Path to functions
addpath(fullfile(projectFolder, 'functions'));


%% =======================================================================
% SETTINGS
%% =======================================================================

nSubjects     = 20;
nSamples      = 10;
nSimulations  = 1000;      % model simulations per trial
nPermutations = 10000;   % permutations of the permutation tests



% Kernels: subjects x samples x RT group (1 = long RT, 2 = short RT)
kernel_selected          = nan(nSubjects,nSamples,2);
kernel_nonSelected       = nan(nSubjects,nSamples,2);
kernel_selected_model    = nan(nSubjects,nSamples,2);
kernel_nonSelected_model = nan(nSubjects,nSamples,2);


%% =======================================================================
% MAIN LOOP ACROSS SUBJECTS
%% =======================================================================

for subject = 1:nSubjects

    %% ===================================================================
    % LOAD SUBJECT DATA
    %
    % All_right / All_left : trials x samples contrasts
    % Real_data            : [choice, confidence, RT, generative mean]
    %% ===================================================================

    load(fullfile( ...
        modelFolder, ...
        ['All_data_model_9_sbj_' num2str(subject) '.mat']));

    nTrials = size(All_right,1);

    realChoice = Real_data(:,1);
    realRT     = Real_data(:,3);

    % Contrast fluctuations around the generative mean of each side
    generativeMean   = Real_data(:,4);
    fluctuationRight = All_right - (0.5 + generativeMean ./ 2);
    fluctuationLeft  = All_left  - (0.5 - generativeMean ./ 2);


    %% ===================================================================
    % MODEL SIMULATION (best-fitting parameters of this subject)
    %
    % modelChoice       : trials x simulations (0 = left, 1 = right)
    % modelDecisionTime : trials x simulations, first bound crossing
    %                     (NaN = bound not reached; these simulations are
    %                     excluded from the median and from both groups)
    %% ===================================================================

    load(fullfile(modelFolder, ['Parameters_model_9_sbj_' num2str(subject) '.mat']));
    bestParameters = Parameters(subject,:,2);

    simulationInput.right = All_right';
    simulationInput.left  = All_left';

    [~, Choices, ~, ~, ~, ~, ~, ~, ~, ~, Cross_time] = ...
        leaky_integrator_simulations_fine( ...
            bestParameters, ...
            simulationInput, ...
            4, ...
            nSimulations);

    modelChoice       = reshape(Choices,    nTrials, nSimulations);
    modelDecisionTime = reshape(Cross_time, nTrials, nSimulations);


    %% ===================================================================
    % RT GROUPS (median split)
    %
    % rtGroup = 1: long RT  (> median)
    % rtGroup = 2: short RT (< median)
    %% ===================================================================

    realRTgroup = nan(nTrials,1);
    realRTgroup(realRT > nanmedian(realRT)) = 1;
    realRTgroup(realRT < nanmedian(realRT)) = 2;

    modelRTgroup = nan(nTrials,nSimulations);
    modelRTgroup(modelDecisionTime > nanmedian(modelDecisionTime(:))) = 1;
    modelRTgroup(modelDecisionTime < nanmedian(modelDecisionTime(:))) = 2;


    %% ===================================================================
    % KERNELS FOR LONG (1) AND SHORT (2) RT
    %
    % Selected side    : right stimulus for right choices,
    %                    left stimulus for left choices
    % Non-selected side: the other stimulus
    % Kernels of right and left choices are averaged.
    %
    % Model: average across all selected (trial, simulation) pairs, i.e.
    % average across trials weighted by the number of selected simulations
    %% ===================================================================

    for rtGroup = 1:2

        % Data
        rightTrials = realChoice == 1 & realRTgroup == rtGroup;
        leftTrials  = realChoice == 0 & realRTgroup == rtGroup;

        kernel_selected(subject,:,rtGroup) = mean([ ...
            mean(fluctuationRight(rightTrials,:),1); ...
            mean(fluctuationLeft(leftTrials,:),1)], 1);

        kernel_nonSelected(subject,:,rtGroup) = mean([ ...
            mean(fluctuationLeft(rightTrials,:),1); ...
            mean(fluctuationRight(leftTrials,:),1)], 1);

        % Model: number of simulations of each trial in this group
        rightCounts = sum(modelChoice == 1 & modelRTgroup == rtGroup, 2);
        leftCounts  = sum(modelChoice == 0 & modelRTgroup == rtGroup, 2);

        kernel_selected_model(subject,:,rtGroup) = mean([ ...
            (rightCounts' * fluctuationRight) ./ sum(rightCounts); ...
            (leftCounts'  * fluctuationLeft)  ./ sum(leftCounts)], 1);

        kernel_nonSelected_model(subject,:,rtGroup) = mean([ ...
            (rightCounts' * fluctuationLeft)  ./ sum(rightCounts); ...
            (leftCounts'  * fluctuationRight) ./ sum(leftCounts)], 1);

    end

end


%% =======================================================================
% KERNEL ASYMMETRY (Figure S3F)
%

%% =======================================================================

asymmetry = ...
    abs(kernel_selected(:,:,2)    - kernel_selected(:,:,1)) - ...
    abs(kernel_nonSelected(:,:,2) - kernel_nonSelected(:,:,1));

asymmetry_model = ...
    abs(kernel_selected_model(:,:,2)    - kernel_selected_model(:,:,1)) - ...
    abs(kernel_nonSelected_model(:,:,2) - kernel_nonSelected_model(:,:,1));

pValues = nan(1,nSamples);

for sample = 1:nSamples
    pValues(sample) = permutationTest( ...
        asymmetry(:,sample), ...
        zeros(nSubjects,1), ...
        nPermutations);
end

significantSamples = double(fdr_bh(pValues,0.05));
significantSamples(significantSamples == 0) = nan;


%% =======================================================================
% FIGURE 1 (Figure S3F) - KERNEL ASYMMETRY
%% =======================================================================

gray = [0.5 0.5 0.5];

figure(1)
hold on

% Data: group average +/- SEM
e = errorbar( ...
    1:nSamples, ...
    nanmean(asymmetry,1), ...
    nanstd(asymmetry,0,1) ./ sqrt(nSubjects), ...
    'color', 'k', ...
    'linewidth', 5, ...
    'linestyle', 'none');
e.CapSize = 0;

scatter( ...
    1:nSamples, ...
    nanmean(asymmetry,1), ...
    80, ...
    'k', ...
    'filled', ...
    'MarkerEdgeColor', [1 1 1]);

% Model: group average +/- SEM
shadedErrorBar( ...
    1:nSamples, ...
    nanmean(asymmetry_model,1), ...
    nanstd(asymmetry_model,0,1) ./ sqrt(nSubjects), ...
    'lineProps', {'.-','color','k'});

% Significant samples (p < 0.05, FDR-corrected)
scatter( ...
    1:nSamples, ...
    significantSamples .* (0.05/3), ...
    80, ...
    'k', ...
    'square', ...
    'filled');

yline(0,'--','color',gray);

xticks(1:nSamples)
xlim([0 nSamples+1])
ylabel({'Confidence kernel asymmetry'; '|1-2| - |4-3|'})
xlabel('Sample number')
set(gca,'fontsize',15)


%% =======================================================================
% FIGURE 2 (Figure S3E) - THE FOUR KERNELS
%
%   1 = selected, long RT          (black)
%   2 = selected, short RT         (gray)
%   3 = non-selected, long RT      (dark blue)
%   4 = non-selected, short RT     (light blue)
%% =======================================================================

% Colors: rows = long / short RT
colors_selected    = [0 0 0; gray];
colors_nonSelected = [0 0.4470 0.7410; [0 162 255] ./ 255];

figure(2)
hold on

for side = 1:2

    if side == 1
        kernel_data  = kernel_selected;
        kernel_model = kernel_selected_model;
        colors       = colors_selected;
    else
        kernel_data  = kernel_nonSelected;
        kernel_model = kernel_nonSelected_model;
        colors       = colors_nonSelected;
    end

    for rtGroup = 1:2

        % Data: group average +/- SEM
        e = errorbar( ...
            1:nSamples, ...
            nanmean(kernel_data(:,:,rtGroup),1), ...
            nanstd(kernel_data(:,:,rtGroup),0,1) ./ sqrt(nSubjects), ...
            'color', colors(rtGroup,:), ...
            'linewidth', 5, ...
            'linestyle', 'none');
        e.CapSize = 0;

        scatter( ...
            1:nSamples, ...
            nanmean(kernel_data(:,:,rtGroup),1), ...
            80, ...
            colors(rtGroup,:), ...
            'filled', ...
            'MarkerEdgeColor', [1 1 1]);

        % Model: group average +/- SEM
        shadedErrorBar( ...
            1:nSamples, ...
            nanmean(kernel_model(:,:,rtGroup),1), ...
            nanstd(kernel_model(:,:,rtGroup),0,1) ./ sqrt(nSubjects), ...
            'lineProps', {'-','color',colors(rtGroup,:)});

    end

end

yline(0,'--','color',gray);

xticks(1:nSamples)
xlim([0 nSamples+1])
ylabel('Excess contrast')
xlabel('Sample number')
set(gca,'fontsize',15)
