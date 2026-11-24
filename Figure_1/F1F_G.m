%% =======================================================================
% Figure 1F-G
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

% Path to model results (behavior + model predictions of each subject)
modelFolder = fullfile(projectFolder, 'data', 'model', 'results');

% Path to functions
addpath(fullfile(projectFolder, 'functions'));


%% =======================================================================
% SETTINGS
%% =======================================================================

nSubjects     = 20;
nSamples      = 10;
nSimulations  = 150;      % simulated choice/confidence reports per trial
nPermutations = 10000;   % permutations of the permutation tests

% Kernels: subjects x samples x confidence (1 = low, 2 = high)
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
    % All_Conf             : model probabilities of
    %                        [right/high, right/low, left/high, left/low]
    % Real_data            : [choice, confidence, RT, generative mean]
    %% ===================================================================

    load(fullfile( ...
        modelFolder, ...
        ['All_data_model_9_sbj_' num2str(subject) '.mat']));

    % Participant choice (1 = right, 0 = left) and confidence (1 = high, 0 = low)
    realChoice     = Real_data(:,1);
    realConfidence = Real_data(:,2) - 1;

    % Contrast fluctuations around the generative mean of each side
    generativeMean   = Real_data(:,4);
    fluctuationRight = All_right - (0.5 + generativeMean ./ 2);
    fluctuationLeft  = All_left  - (0.5 - generativeMean ./ 2);


    %% ===================================================================
    % SIMULATED CHOICE/CONFIDENCE REPORTS
    %
    % simulatedCounts(trial,k): number of simulations (out of
    % nSimulations) in which the model gave report k
    % k: 1 = right/high, 2 = right/low, 3 = left/high, 4 = left/low
    %% ===================================================================

    modelProbabilities = abs(All_Conf) ./ sum(abs(All_Conf),2);
    simulatedCounts    = mnrnd(nSimulations, modelProbabilities);


    %% ===================================================================
    % KERNELS FOR LOW (1) AND HIGH (2) CONFIDENCE
    %
    
    %% ===================================================================

    for confidenceLevel = 1:2

        % Data
        rightTrials = realChoice == 1 & realConfidence == confidenceLevel - 1;
        leftTrials  = realChoice == 0 & realConfidence == confidenceLevel - 1;

        kernel_selected(subject,:,confidenceLevel) = mean([ ...
            mean(fluctuationRight(rightTrials,:),1); ...
            mean(fluctuationLeft(leftTrials,:),1)], 1);

        kernel_nonSelected(subject,:,confidenceLevel) = mean([ ...
            mean(fluctuationLeft(rightTrials,:),1); ...
            mean(fluctuationRight(leftTrials,:),1)], 1);

        % Model (low confidence: reports 2 and 4, high: reports 1 and 3)
        if confidenceLevel == 1
            rightCounts = simulatedCounts(:,2);
            leftCounts  = simulatedCounts(:,4);
        else
            rightCounts = simulatedCounts(:,1);
            leftCounts  = simulatedCounts(:,3);
        end

        kernel_selected_model(subject,:,confidenceLevel) = mean([ ...
            (rightCounts' * fluctuationRight) ./ sum(rightCounts); ...
            (leftCounts'  * fluctuationLeft)  ./ sum(leftCounts)], 1);

        kernel_nonSelected_model(subject,:,confidenceLevel) = mean([ ...
            (rightCounts' * fluctuationLeft)  ./ sum(rightCounts); ...
            (leftCounts'  * fluctuationRight) ./ sum(leftCounts)], 1);

    end

end


%% =======================================================================
% CONFIDENCE KERNEL ASYMMETRY (Figure 1G)
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
% FIGURE 1 (Figure 1G) - CONFIDENCE KERNEL ASYMMETRY
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
% FIGURE 2 (Figure 1F) - THE FOUR KERNELS
%
%   1 = selected, high confidence        (gray)
%   2 = selected, low confidence         (black)
%   3 = non-selected, low confidence     (dark blue)
%   4 = non-selected, high confidence    (light blue)
%% =======================================================================

% Colors: rows = low / high confidence
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

    for confidenceLevel = 1:2

        % Data: group average +/- SEM
        e = errorbar( ...
            1:nSamples, ...
            nanmean(kernel_data(:,:,confidenceLevel),1), ...
            nanstd(kernel_data(:,:,confidenceLevel),0,1) ./ sqrt(nSubjects), ...
            'color', colors(confidenceLevel,:), ...
            'linewidth', 5, ...
            'linestyle', 'none');
        e.CapSize = 0;

        scatter( ...
            1:nSamples, ...
            nanmean(kernel_data(:,:,confidenceLevel),1), ...
            80, ...
            colors(confidenceLevel,:), ...
            'filled', ...
            'MarkerEdgeColor', [1 1 1]);

        % Model: group average +/- SEM
        shadedErrorBar( ...
            1:nSamples, ...
            nanmean(kernel_model(:,:,confidenceLevel),1), ...
            nanstd(kernel_model(:,:,confidenceLevel),0,1) ./ sqrt(nSubjects), ...
            'lineProps', {'-','color',colors(confidenceLevel,:)});

    end

end

yline(0,'--','color',gray);

xticks(1:nSamples)
xlim([0 nSamples+1])
ylabel('Excess contrast')
xlabel('Sample number')
set(gca,'fontsize',15)
