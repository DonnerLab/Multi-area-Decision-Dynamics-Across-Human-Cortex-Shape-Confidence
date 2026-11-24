%% =======================================================================
% Figure 1D-E
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

% Data folders
behaviorFolder = fullfile(projectFolder, 'data', 'behavior');
modelFolder    = fullfile(projectFolder, 'data', 'model', 'results');

% Path to functions
addpath(fullfile(projectFolder, 'functions'));

% Behavioral data of all subjects (All_Data: subjects x drug conditions)
load(fullfile(behaviorFolder, 'data.mat'))


%% =======================================================================
% SETTINGS
%% =======================================================================

nSubjects     = 20;
nBins         = 4;        % evidence-strength quantile bins
nSimulations  = 250;      % model simulations per trial
nPermutations = 10000;   % permutations of the permutation tests


green = [0.4660 0.6740 0.1880];

% subjects x bins; third dimension: 1 = correct, 2 = error trials
meanConfidence       = nan(nSubjects,nBins,2);
meanEvidence         = nan(nSubjects,nBins,2);
meanConfidence_model = nan(nSubjects,nBins,2);
meanEvidence_model   = nan(nSubjects,nBins,2);

% subjects x 2 (1 = correct, 2 = error trials)
correlation_confidence_evidence = nan(nSubjects,2);


%% =======================================================================
% MAIN LOOP ACROSS SUBJECTS
%% =======================================================================

for subject = 1:nSubjects

    %% ===================================================================
    % LOAD SUBJECT DATA (the 3 drug conditions combined)
    %
    % columns:  1      = choice (0 = left, 1 = right)
    %           2:11   = left contrast samples
    %           12:21  = right contrast samples
    %           22     = confidence (0 = low, 1 = high)
    %% ===================================================================

    subjectData = [];

    for drugCondition = 1:3
        subjectData = [subjectData; cell2mat(All_Data(subject,drugCondition))];
    end

    % Remove trials without a response
    subjectData(isnan(subjectData(:,1)) | isnan(subjectData(:,22)),:) = [];

    choice           = subjectData(:,1);
    leftStimulus     = subjectData(:,2:11);
    rightStimulus    = subjectData(:,12:21);
    confidenceReport = subjectData(:,22);

    % Correct side (1 = right, 0 = left), accuracy and evidence strength
    correctSide      = double(mean(rightStimulus,2) > mean(leftStimulus,2));
    accuracy         = double(choice == correctSide);
    evidenceStrength = abs(mean(rightStimulus,2) - mean(leftStimulus,2));


    %% ===================================================================
    % MODEL SIMULATION (best-fitting parameters of this subject)
    %
  
    %% ===================================================================

    load(fullfile(modelFolder, ['Parameters_model_9_sbj_' num2str(subject) '.mat']));
    bestParameters = Parameters(subject,:,2);

    simulationInput.t     = 10;
    simulationInput.right = rightStimulus';
    simulationInput.left  = leftStimulus';

    [~, ~, ~, modelChoice, ~, modelConfidence] = ...
        model_simulation(bestParameters, simulationInput, 4, nSimulations);

    nTrials = length(choice);

    modelChoice     = reshape(modelChoice,     nTrials, nSimulations);
    modelConfidence = reshape(modelConfidence, nTrials, nSimulations);

    modelAccuracy = double(modelChoice == repmat(correctSide,1,nSimulations));
    modelEvidence = repmat(evidenceStrength,1,nSimulations);

    modelConfidence = modelConfidence(:);
    modelAccuracy   = modelAccuracy(:);
    modelEvidence   = modelEvidence(:);


    %% ===================================================================
    % EVIDENCE BINS
    %
    % Quantiles of the participant's evidence strength, applied to the
    % data and to the model simulations
    %% ===================================================================

    binEdges = [-inf quantile(evidenceStrength, nBins-1) inf];

    [~, ~, evidenceBin]       = histcounts(evidenceStrength, binEdges);
    [~, ~, evidenceBin_model] = histcounts(modelEvidence,    binEdges);


    %% ===================================================================
    % CORRECT (outcome = 1) AND ERROR (outcome = 2) TRIALS
    %% ===================================================================

    for outcome = 1:2

        isOutcome       = accuracy      == (outcome == 1);
        isOutcome_model = modelAccuracy == (outcome == 1);

        % Correlation between confidence and evidence strength (Figure 1E)
        correlation_confidence_evidence(subject,outcome) = corr( ...
            confidenceReport(isOutcome), ...
            evidenceStrength(isOutcome), ...
            'type', 'spearman');

        % Average per evidence bin (Figure 1D)
        for bin = 1:nBins

            meanConfidence(subject,bin,outcome) = ...
                nanmean(confidenceReport(isOutcome & evidenceBin == bin));
            meanEvidence(subject,bin,outcome) = ...
                nanmean(evidenceStrength(isOutcome & evidenceBin == bin));

            meanConfidence_model(subject,bin,outcome) = ...
                nanmean(modelConfidence(isOutcome_model & evidenceBin_model == bin));
            meanEvidence_model(subject,bin,outcome) = ...
                nanmean(modelEvidence(isOutcome_model & evidenceBin_model == bin));

        end

    end

end


%% =======================================================================
% STATISTICS (Figure 1E)

%% =======================================================================

pValue_correct = permutationTest( ...
    correlation_confidence_evidence(:,1), zeros(nSubjects,1), nPermutations);

pValue_error = permutationTest( ...
    correlation_confidence_evidence(:,2), zeros(nSubjects,1), nPermutations);

pValue_correct_vs_error = permutationTest( ...
    correlation_confidence_evidence(:,1) - correlation_confidence_evidence(:,2), ...
    zeros(nSubjects,1), nPermutations);


%% =======================================================================
% FIGURE 1 (Figure 1D) - CONFIDENCE VS EVIDENCE, CORRECT / ERROR TRIALS
%% =======================================================================

outcomeColors = [green; 1 0 0];

figure(1)
hold on

for outcome = 1:2

    % Model: group average +/- SEM
    shadedErrorBar( ...
        nanmean(meanEvidence_model(:,:,outcome),1), ...
        nanmean(meanConfidence_model(:,:,outcome),1) .* 100, ...
        nanstd(meanConfidence_model(:,:,outcome) .* 100,0,1) ./ sqrt(nSubjects), ...
        'lineProps', {'-','color',outcomeColors(outcome,:)});

    % Data: group average
    scatter( ...
        nanmean(meanEvidence(:,:,outcome),1), ...
        nanmean(meanConfidence(:,:,outcome),1) .* 100, ...
        70, ...
        outcomeColors(outcome,:), ...
        'filled', ...
        'MarkerEdgeColor', [1 1 1]);

end

xlabel('| \Delta Contrast (R-L) |')
ylabel('High confidence report (%)')
set(gca,'fontsize',15)


%% =======================================================================
% FIGURE 2 (Figure 1E) - CORRELATION OF CONFIDENCE WITH EVIDENCE
%% =======================================================================

figure(2)
hold on

for outcome = 1:2

    bar( ...
        outcome, ...
        mean(correlation_confidence_evidence(:,outcome)), ...
        'facecolor', outcomeColors(outcome,:), ...
        'facealpha', 0.5);

    for subject = 1:nSubjects
        scatter( ...
            outcome + (2*rand(1) - 1) * 0.05, ...
            correlation_confidence_evidence(subject,outcome), ...
            100, ...
            outcomeColors(outcome,:), ...
            'filled', ...
            'markerfacealpha', 0.5);
    end

end

yline(0, 'linewidth', 3, 'color', [0.5 0.5 0.5], 'linestyle', '--')

xticks([])
ylabel({'Scaling of confidence report'; 'with evidence'})
title({sprintf('correct: p = %.2g,  error: p = %.2g', pValue_correct, pValue_error); ...
       sprintf('correct vs error: p = %.2g', pValue_correct_vs_error)})
set(gca,'fontsize',15)
