%% =======================================================================
% Figure 1B and Figure S3A-B
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
nBins         = 5;        % evidence-strength quantile bins
nSimulations  = 1000;      % model simulations per trial
nPermutations = 10000;   % permutations of the permutation test

% Duration of one stimulus sample (s): converts the model crossing time
% from units of samples to seconds
sampleDuration = 0.1;

% Simulations in which the bound was not reached get the crossing time of
% the end of the stimulus (10 samples)
noCrossingTime = 10;

meanEvidence           = nan(nSubjects,nBins);
meanConfidence         = nan(nSubjects,nBins);
meanConfidence_model   = nan(nSubjects,nBins);
meanRT                 = nan(nSubjects,nBins);
meanDecisionTime_model = nan(nSubjects,nBins);

correlation_RT_confidence = nan(nSubjects,1);


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
    %           24     = reaction time (s)
    %% ===================================================================

    subjectData = [];

    for drugCondition = 1:3
        subjectData = [subjectData; cell2mat(All_Data(subject,drugCondition))];
    end

    leftStimulus     = subjectData(:,2:11);
    rightStimulus    = subjectData(:,12:21);
    confidenceReport = subjectData(:,22);
    reactionTime     = subjectData(:,24);

    % Evidence strength = |mean right contrast - mean left contrast|
    evidenceStrength = abs(mean(rightStimulus,2) - mean(leftStimulus,2));


    %% ===================================================================
    % MODEL SIMULATION (best-fitting parameters of this subject)
    %
    % Confidence : trials x 4 probabilities
    %              [right/high, right/low, left/high, left/low]
    % Cross_time : 1 x trials x simulations, first bound crossing
    %              (in samples; NaN = bound not reached)
    %% ===================================================================

    load(fullfile(modelFolder, ['Parameters_model_9_sbj_' num2str(subject) '.mat']));
    bestParameters = Parameters(subject,:,2);

    simulationInput.right = rightStimulus';
    simulationInput.left  = leftStimulus';

    [~, ~, Confidence, ~, ~, ~, ~, ~, ~, ~, Cross_time] = ...
        leaky_integrator_simulations_fine( ...
            bestParameters, ...
            simulationInput, ...
            4, ...
            nSimulations);

    Cross_time(isnan(Cross_time)) = noCrossingTime;

    % Model probability of a high-confidence report
    highConfidence_model = Confidence(:,1) + Confidence(:,3);

    % Model decision time (s): crossing time averaged across simulations
    decisionTime_model = squeeze(mean(Cross_time,3))' .* sampleDuration;


    %% ===================================================================
    % CORRELATION BETWEEN RT AND CONFIDENCE (Figure S3A)
    %% ===================================================================

    correlation_RT_confidence(subject) = corr( ...
        reactionTime, ...
        confidenceReport, ...
        'type', 'spearman', ...
        'rows', 'complete');


    %% ===================================================================
    % AVERAGE PER EVIDENCE-STRENGTH BIN (5 quantiles)
    %% ===================================================================

    evidenceQuantiles = quantile(evidenceStrength, nBins-1);

    [~, ~, evidenceBin] = histcounts( ...
        evidenceStrength, ...
        [-inf evidenceQuantiles inf]);

    for bin = 1:nBins

        meanEvidence(subject,bin)           = nanmean(evidenceStrength(evidenceBin == bin));
        meanConfidence(subject,bin)         = nanmean(confidenceReport(evidenceBin == bin));
        meanConfidence_model(subject,bin)   = nanmean(highConfidence_model(evidenceBin == bin));
        meanRT(subject,bin)                 = nanmean(reactionTime(evidenceBin == bin));
        meanDecisionTime_model(subject,bin) = nanmean(decisionTime_model(evidenceBin == bin));

    end

end


%% =======================================================================
% FIGURE 1 (Figure 1B) - HIGH-CONFIDENCE REPORTS VS EVIDENCE STRENGTH
%% =======================================================================

gray = [0.5 0.5 0.5];

figure(1)
hold on

% Data
e = errorbar( ...
    mean(meanEvidence,1), ...
    mean(meanConfidence,1) .* 100, ...
    std(meanConfidence .* 100,0,1) ./ sqrt(nSubjects), ...
    'color', gray, ...
    'linewidth', 5, ...
    'linestyle', 'none');
e.CapSize = 0;

scatter( ...
    mean(meanEvidence,1), ...
    mean(meanConfidence,1) .* 100, ...
    70, ...
    gray, ...
    'filled', ...
    'MarkerEdgeColor', [1 1 1]);

% Model
shadedErrorBar( ...
    mean(meanEvidence,1), ...
    mean(meanConfidence_model,1) .* 100, ...
    std(meanConfidence_model .* 100,0,1) ./ sqrt(nSubjects), ...
    'lineProps', {'-','color',gray});

xlabel('| \Delta Contrast (R-L) |')
ylabel('High confidence report (%)')
set(gca,'fontsize',15)


%% =======================================================================
% FIGURE 2 (Figure S3A) - CORRELATION BETWEEN CONFIDENCE AND RT
%
%% =======================================================================

pValue_RT_confidence = permutationTest( ...
    correlation_RT_confidence, ...
    zeros(nSubjects,1), ...
    nPermutations);

figure(2)
hold on

bar(1, mean(correlation_RT_confidence), 'facecolor', gray);

for subject = 1:nSubjects
    scatter( ...
        1 + (rand(1)*2 - 1) * 0.05, ...
        correlation_RT_confidence(subject), ...
        100, ...
        'k', ...
        'filled', ...
        'markerfacealpha', 0.5);
end

yline(0, 'linewidth', 3, 'color', gray, 'linestyle', '--')

xticks([])
ylabel({'Correlation coefficient'; 'confidence with RT'})
title(sprintf('p = %.2g', pValue_RT_confidence))
set(gca,'fontsize',15)


%% =======================================================================
% FIGURE 3 (Figure S3B) - RT VS EVIDENCE STRENGTH
%
% Left y-axis : RT (data)
% Right y-axis: model decision time
%% =======================================================================

figure(3)

yyaxis left
hold on

e = errorbar( ...
    mean(meanEvidence,1), ...
    nanmean(meanRT,1), ...
    nanstd(meanRT,0,1) ./ sqrt(nSubjects), ...
    'color', 'k', ...
    'linewidth', 5, ...
    'linestyle', 'none');
e.CapSize = 0;

scatter( ...
    mean(meanEvidence,1), ...
    nanmean(meanRT,1), ...
    70, ...
    'k', ...
    'filled', ...
    'MarkerEdgeColor', [1 1 1]);

ylabel('Reaction time (s)')
set(gca,'ycolor','k')

yyaxis right
hold on

shadedErrorBar( ...
    mean(meanEvidence,1), ...
    nanmean(meanDecisionTime_model,1), ...
    nanstd(meanDecisionTime_model,0,1) ./ sqrt(nSubjects), ...
    'lineProps', {'-','color','k'});

ylabel('Model decision time (s)')
set(gca,'ycolor','k')

xlabel('| \Delta Contrast (R-L) |')
set(gca,'fontsize',15)
