%% =======================================================================
% Plot Figure 1B-C and Figure S1 G
%% =======================================================================

clear all
close all

% Folder where this script is located
scriptFolder = fileparts(mfilename('fullpath'));

% Go up two folders
projectFolder = fileparts(fileparts(scriptFolder));

% Path to behavior data
dataFolder = fullfile(projectFolder, 'data', 'behavior');

% Path to functions
functionsFolder = fullfile(projectFolder, 'functions');
addpath(functionsFolder)

% Full path to data file
filePath = fullfile(dataFolder, 'data.mat');

% Load file
load(filePath)

% Model parameters
modelFolder = fullfile(projectFolder, 'data', 'model', 'results');

% Path to model_simulation.m
functionsFolder = fullfile(projectFolder, 'functions');
addpath(simulationFolder);


%% =======================================================================
% RT and High Confidence vs Evidence Strength
% for model predictions and behavioral data
%% =======================================================================

for subject = 1:20
    
    %% ===================================================================
    % LOAD SUBJECT DATA
    %% ===================================================================
    
    subjectData = [];
    
    % Combine the 3 drug conditions
    for drugCondition = 1:3
        
        currentBlock = cell2mat(All_Data(subject, drugCondition));
        subjectData  = [subjectData; currentBlock];
        
    end
    
    
    %% ===================================================================
    % EXTRACT BEHAVIORAL VARIABLES
    %% ===================================================================
    
    % Subject choices
    responseChoice = subjectData(:,1);
    
    % Stimulus contrasts
    stimulusMatrix = subjectData(:,2:21);
    
    % Reaction time
    reactionTime = subjectData(:,24);
    
    % Confidence reports
    confidenceReport = subjectData(:,22);
    
    
    %% ===================================================================
    % COMPUTE STIMULUS FEATURES
    %% ===================================================================
    
    % Mean evidence on each side
    meanRightStimulus = mean(stimulusMatrix(:,1:10), 2);
    meanLeftStimulus  = mean(stimulusMatrix(:,11:20), 2);
    
    % True stimulus side
    trueStimulusSide = double((meanRightStimulus - meanLeftStimulus) < 0);
    trueStimulusSide(trueStimulusSide == 0) = 0;
    
    % Evidence strength
    evidenceStrength = abs(meanRightStimulus - meanLeftStimulus);
    
    
    %% ===================================================================
    % LOAD BEST MODEL PARAMETERS
    %% ===================================================================
    
    clearvars Parameters bestParameters
    
    load(fullfile( ...
        modelFolder, ...
        ['Parameters_model_9_sbj_' num2str(subject) '.mat']));
    
    bestParameters = Parameters(subject,:,2);
    
    
    %% ===================================================================
    % MODEL SIMULATION
    %% ===================================================================
    
    simulationInput.t     = 10;
    simulationInput.right = stimulusMatrix(:,1:10)';
    simulationInput.left  = stimulusMatrix(:,11:20)';
    simulationInput.rt    = reactionTime;
    
    [decisionVariable, probabilityLeft, probabilityRight, ...
     simulatedChoices, simulatedConfidence, modelOutput, ...
     collapsingThreshold, confidenceAtDecision, ...
     winningConfidence, losingConfidence] = ...
        model_simulation(bestParameters, simulationInput, 4);
    
    
    %% ===================================================================
    % MODEL OUTPUTS
    %% ===================================================================
    
    % Final model confidence per trial
    modelConfidence = confidenceAtDecision;
    
    
    %% ===================================================================
    % EVIDENCE BINS (5 QUANTILES)
    %% ===================================================================
    
    evidenceQuantiles = quantile(evidenceStrength(:,1), 4);
    
    [~, ~, evidenceBin] = histcounts( ...
        evidenceStrength(:,1), ...
        [-inf evidenceQuantiles inf]);
    
    
    %% ===================================================================
    % RT PROXY FROM MODEL CONFIDENCE
    %% ===================================================================
    
    % Real RT
    realRT = reactionTime;
    
    % Transform model confidence into RT proxy
    modelRT = 1 - (mean(modelConfidence,3))';
    
    % Align model RT mean with empirical RT mean
    modelRT = modelRT + ...
        (nanmean(realRT) - nanmean(modelRT));
    
    
    %% ===================================================================
    % HIGH CONFIDENCE REPORTS
    %% ===================================================================
    
    % Real confidence report
    realConfidence = confidenceReport;
    
    % Model confidence report
    Conf_model = sum([simulatedConfidence(:,1), ...
                      simulatedConfidence(:,3)], 2);
    
    
    %% ===================================================================
    % CORRELATION BETWEEN RT AND CONFIDENCE
    %% ===================================================================
    
    [pValue, correlationValue] = corr( ...
        realRT, ...
        confidenceReport, ...
        'type', 'spearman', ...
        'rows', 'complete');
    
    correlation_RT_confidence(subject) = pValue;
    
    
    %% ===================================================================
    % COMPUTE MEAN VALUES PER EVIDENCE BIN
    %% ===================================================================
    
    for bin = 1:5
        
        % Real data
        meanRT_real(subject,bin) = ...
            nanmean(realRT(evidenceBin == bin,1));
        
        meanConfidence_real(subject,bin) = ...
            nanmean(realConfidence(evidenceBin == bin,1));
        
        meanEvidence_real(subject,bin) = ...
            nanmean(evidenceStrength(evidenceBin == bin,1));
        
        
        % Model
        meanRT_model(subject,bin) = ...
            nanmean(modelRT(evidenceBin == bin,1));
        
        meanConfidence_model(subject,bin) = ...
            nanmean(Conf_model(evidenceBin == bin,1));
        
        meanEvidence_model(subject,bin) = ...
            nanmean(evidenceStrength(evidenceBin == bin,1));
        
    end
    
end


%% =======================================================================
% RESHAPE FOR PLOTTING
%% =======================================================================

RT_real   = reshape(meanRT_real,         20, 5);
EV_real   = reshape(meanEvidence_real,   20, 5);
CONF_real = reshape(meanConfidence_real, 20, 5) .* 100;

RT_model   = reshape(meanRT_model,         20, 5);
EV_model   = reshape(meanEvidence_model,   20, 5);
CONF_model = reshape(meanConfidence_model, 20, 5) .* 100;


%% =======================================================================
% FIGURE 1 - REACTION TIME
%% =======================================================================

figure(1)
clf
hold on

% Real data
e = errorbar( ...
    mean(EV_real,1), ...
    mean(RT_real,1), ...
    std(RT_real)./sqrt(20), ...
    'color', 'k', ...
    'linewidth', 5, ...
    'linestyle', 'none');

e.CapSize = 0;

scatter( ...
    mean(EV_real,1), ...
    mean(RT_real,1), ...
    70, ...
    'k', ...
    'filled', ...
    'MarkerEdgeColor', [1 1 1]);

% Model
shadedErrorBar( ...
    mean(EV_real,1), ...
    mean(RT_model,1), ...
    std(RT_model)./sqrt(20), ...
    'lineProps', {'-','color','k'});

xlabel('Evidence strength')
ylabel('Reaction Time (s)')
set(gca,'fontsize',15)


%% =======================================================================
% FIGURE 2 - HIGH CONFIDENCE REPORTS
%% =======================================================================

figure(2)
clf
hold on

% Real data
e = errorbar( ...
    mean(EV_real,1), ...
    mean(CONF_real,1), ...
    std(CONF_real)./sqrt(20), ...
    'color', [0.5 0.5 0.5], ...
    'linewidth', 5, ...
    'linestyle', 'none');

e.CapSize = 0;

scatter( ...
    mean(EV_real,1), ...
    mean(CONF_real,1), ...
    70, ...
    [0.5 0.5 0.5], ...
    'filled', ...
    'MarkerEdgeColor', [1 1 1]);

% Model
shadedErrorBar( ...
    mean(EV_real,1), ...
    mean(CONF_model,1), ...
    std(CONF_model)./sqrt(20), ...
    'lineProps', {'-','color',[0.5 0.5 0.5]});

xlabel('Evidence strength')
ylabel('Prob. High Confidence')
set(gca,'fontsize',15)


%% =======================================================================
% FIGURE 3 - CORRELATION BETWEEN RT AND CONFIDENCE
%% =======================================================================

figure(3)
clf
hold on

bar( ...
    1, ...
    mean(correlation_RT_confidence,2), ...
    'facecolor', [0.5 0.5 0.5]);

for subject = 1:20
    
    scatter( ...
        1 + ((rand(1)*2 - 1) * 0.05), ...
        correlation_RT_confidence(subject), ...
        100, ...
        'k', ...
        'filled', ...
        'markerfacealpha', 0.5);
    
end

xticks([])

yline(0, ...
    'linewidth', 3, ...
    'color', [0.5 0.5 0.5], ...
    'linestyle', '--')

ylabel({'Correlation Coefficient'; ...
         'Confidence with RT'})

set(gca,'fontsize',15)