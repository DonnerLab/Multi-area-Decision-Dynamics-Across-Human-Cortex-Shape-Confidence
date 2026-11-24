%% =======================================================================
% Figure S3C-D
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
nSimulations  = 1000;      % model simulations per trial
nPermutations = 10000;   % permutations of the permutation tests

% Duration of one stimulus sample (s)
sampleDuration = 0.1;


green = [0.4660 0.6740 0.1880];

% subjects x bins; third dimension: 1 = correct, 2 = error trials
meanRT                 = nan(nSubjects,nBins,2);
meanEvidence           = nan(nSubjects,nBins,2);
meanDecisionTime_model = nan(nSubjects,nBins,2);
meanEvidence_model     = nan(nSubjects,nBins,2);

% subjects x 2 (1 = correct, 2 = error trials)
correlation_RT_evidence = nan(nSubjects,2);


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
    %           24     = reaction time (s)
    %% ===================================================================

    subjectData = [];

    for drugCondition = 1:3
        subjectData = [subjectData; cell2mat(All_Data(subject,drugCondition))];
    end

    % Remove trials without a response
    subjectData(isnan(subjectData(:,1)) | isnan(subjectData(:,24)),:) = [];

    choice        = subjectData(:,1);
    leftStimulus  = subjectData(:,2:11);
    rightStimulus = subjectData(:,12:21);
    reactionTime  = subjectData(:,24);

    % Correct side (1 = right, 0 = left), accuracy and evidence strength
    correctSide      = double(mean(rightStimulus,2) > mean(leftStimulus,2));
    accuracy         = double(choice == correctSide);
    evidenceStrength = abs(mean(rightStimulus,2) - mean(leftStimulus,2));


    %% ===================================================================
    % MODEL SIMULATION (best-fitting parameters of this subject)
    %
    % Every simulation of every trial is one model trial
    % (trials x simulations, vectorized)
    % Cross_time: first bound crossing (in samples; NaN = not reached)
    %% ===================================================================

    load(fullfile(modelFolder, ['Parameters_model_9_sbj_' num2str(subject) '.mat']));
    bestParameters = Parameters(subject,:,2);

    simulationInput.right = rightStimulus';
    simulationInput.left  = leftStimulus';

    [~, Choices, ~, ~, ~, ~, ~, ~, ~, ~, Cross_time] = ...
        leaky_integrator_simulations_fine( ...
            bestParameters, ...
            simulationInput, ...
            4, ...
            nSimulations);

    nTrials = length(choice);

    modelChoice       = reshape(Choices,    nTrials, nSimulations);
    modelDecisionTime = reshape(Cross_time, nTrials, nSimulations) .* sampleDuration;

    modelAccuracy = double(modelChoice == repmat(correctSide,1,nSimulations));
    modelEvidence = repmat(evidenceStrength,1,nSimulations);

    modelDecisionTime = modelDecisionTime(:);
    modelAccuracy     = modelAccuracy(:);
    modelEvidence     = modelEvidence(:);


    %% ===================================================================
    % EVIDENCE BINS
    %
  
    %% ===================================================================

    binEdges = [-inf quantile(evidenceStrength, nBins-1) inf];

    [~, ~, evidenceBin]       = histcounts(evidenceStrength, binEdges);
    [~, ~, evidenceBin_model] = histcounts(modelEvidence,    binEdges);


    %% ===================================================================
    % CORRECT (outcome = 1) AND ERROR (outcome = 2) TRIALS
    %% ===================================================================

    for outcome = 1:2

        isOutcome       = accuracy == (outcome == 1);
        isOutcome_model = modelAccuracy == (outcome == 1) & ~isnan(modelDecisionTime);

        % Correlation between RT and evidence strength (Figure S3D)
        correlation_RT_evidence(subject,outcome) = corr( ...
            reactionTime(isOutcome), ...
            evidenceStrength(isOutcome), ...
            'type', 'spearman');

        % Average per evidence bin (Figure S3C)
        for bin = 1:nBins

            meanRT(subject,bin,outcome) = ...
                nanmean(reactionTime(isOutcome & evidenceBin == bin));
            meanEvidence(subject,bin,outcome) = ...
                nanmean(evidenceStrength(isOutcome & evidenceBin == bin));

            meanDecisionTime_model(subject,bin,outcome) = ...
                nanmean(modelDecisionTime(isOutcome_model & evidenceBin_model == bin));
            meanEvidence_model(subject,bin,outcome) = ...
                nanmean(modelEvidence(isOutcome_model & evidenceBin_model == bin));

        end

    end

end


%% =======================================================================
% STATISTICS 
%% =======================================================================

pValue_correct = permutationTest( ...
    correlation_RT_evidence(:,1), zeros(nSubjects,1), nPermutations);

pValue_error = permutationTest( ...
    correlation_RT_evidence(:,2), zeros(nSubjects,1), nPermutations);

pValue_correct_vs_error = permutationTest( ...
    correlation_RT_evidence(:,1) - correlation_RT_evidence(:,2), ...
    zeros(nSubjects,1), nPermutations);


%% =======================================================================
% FIGURE 1 (Figure S3C) - RT VS EVIDENCE, CORRECT / ERROR TRIALS
%
% Left y-axis : RT (data)
% Right y-axis: model decision time
%% =======================================================================

outcomeColors = [green; 1 0 0];

figure(1)

yyaxis left
hold on

for outcome = 1:2
    scatter( ...
        nanmean(meanEvidence(:,:,outcome),1), ...
        nanmean(meanRT(:,:,outcome),1), ...
        70, ...
        outcomeColors(outcome,:), ...
        'filled', ...
        'MarkerEdgeColor', [1 1 1]);
end

ylabel('Reaction time (s)')
set(gca,'ycolor','k')

yyaxis right
hold on

for outcome = 1:2
    shadedErrorBar( ...
        nanmean(meanEvidence_model(:,:,outcome),1), ...
        nanmean(meanDecisionTime_model(:,:,outcome),1), ...
        nanstd(meanDecisionTime_model(:,:,outcome),0,1) ./ sqrt(nSubjects), ...
        'lineProps', {'-','color',outcomeColors(outcome,:)});
end

ylabel('Model decision time (s)')
set(gca,'ycolor','k')

xlabel('| \Delta Contrast (R-L) |')
set(gca,'fontsize',15)


%% =======================================================================
% FIGURE 2 (Figure S3D) - CORRELATION OF RT WITH EVIDENCE
%% =======================================================================

figure(2)
hold on

for outcome = 1:2

    bar( ...
        outcome, ...
        nanmean(correlation_RT_evidence(:,outcome)), ...
        'facecolor', outcomeColors(outcome,:), ...
        'facealpha', 0.5);

    for subject = 1:nSubjects
        scatter( ...
            outcome + (2*rand(1) - 1) * 0.05, ...
            correlation_RT_evidence(subject,outcome), ...
            100, ...
            outcomeColors(outcome,:), ...
            'filled', ...
            'markerfacealpha', 0.5);
    end

end

yline(0, 'linewidth', 3, 'color', [0.5 0.5 0.5], 'linestyle', '--')

xticks([])
ylabel('Scaling of RT with evidence')
title({sprintf('correct: p = %.2g,  error: p = %.2g', pValue_correct, pValue_error); ...
       sprintf('correct vs error: p = %.2g', pValue_correct_vs_error)})
set(gca,'fontsize',15)
