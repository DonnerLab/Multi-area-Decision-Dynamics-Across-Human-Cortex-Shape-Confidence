%% =======================================================================
% Plot Figure 1D-E and Figure S1 H-I

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


%% ========================================================================
% ANALYSIS LOOP
%
% measureType = 1 --> Reaction Time (RT)
% measureType = 2 --> Confidence
%% ========================================================================

for measureType = 1:2
    
    for subject = 1:20
        
        %% =================================================================
        % LOAD SUBJECT DATA
        %% =================================================================
        
        subjectData = [];
        
        % Combine the 3 drug conditions for each subject
        for drugCondition = 1:3
            
            currentBlock = cell2mat(All_Data(subject, drugCondition));
            subjectData  = [subjectData; currentBlock];
            
        end
        
        
        %% =================================================================
        % EXTRACT MAIN VARIABLES
        %% =================================================================
        
        responseChoice = subjectData(:,1);
        stimulusMatrix = subjectData(:,2:21);
        
        % measureType:
        % 1 = RT
        % 2 = Confidence
        if measureType == 1
            behavioralMeasure = subjectData(:,24); % RT
        else
            behavioralMeasure = subjectData(:,22); % Confidence
        end
        
        
        %% =================================================================
        % COMPUTE STIMULUS FEATURES
        %% =================================================================
        
        % Mean contrast on each side
        meanRightStimulus = mean(stimulusMatrix(:,1:10), 2);
        meanLeftStimulus  = mean(stimulusMatrix(:,11:20), 2);
        
        % Stimulus side
        % 0 = left stronger
        % 1 = right stronger
        trueStimulusSide = double((meanRightStimulus - meanLeftStimulus) < 0);
        trueStimulusSide(trueStimulusSide == 0) = 0;
        
        % Correct vs incorrect choice
        trialAccuracy = double(trueStimulusSide == responseChoice);
        
        % Evidence strength
        evidenceStrength = abs(meanRightStimulus - meanLeftStimulus);
        
        
        %% =================================================================
        % LOAD MODEL PARAMETERS
        %% =================================================================
        
        clearvars Parameters bestParameters
        
        load(fullfile( ...
            modelFolder, ...
            ['Parameters_model_9_sbj_' num2str(subject) '.mat']));
        
        % Best-fitting parameters for this subject
        bestParameters = Parameters(subject,:,2);
        
        
        %% =================================================================
        % RUN MODEL SIMULATION
        %% =================================================================
        
        simulationInput.t = 10;
        simulationInput.right = stimulusMatrix(:,1:10)';
        simulationInput.left  = stimulusMatrix(:,11:20)';
        
        [decisionVariable, probabilityLeft, probabilityRight, ...
         simulatedChoices, simulatedConfidence, modelOutput, ...
         collapsingThreshold, confidenceAtDecision, ...
         winningConfidence, losingConfidence] = ...
            model_simulation(bestParameters, simulationInput, 4);
        
        
        %% =================================================================
        % PREPARE MODEL VARIABLES
        %% =================================================================
        
        modelEvidence = repmat(evidenceStrength, [1,1,250]);
        modelEvidence = permute(modelEvidence, [2,1,3]);
        
        modelStimulusSide = repmat(trueStimulusSide, [1,1,250]);
        modelStimulusSide = permute(modelStimulusSide, [2,1,3]);
        
        
        %% =================================================================
        % REMOVE NaN TRIALS
        %% =================================================================
        
        nanTrials = isnan(responseChoice);
        
        responseChoice(nanTrials)     = [];
        evidenceStrength(nanTrials,:) = [];
        behavioralMeasure(nanTrials)  = [];
        trialAccuracy(nanTrials)      = [];
        trueStimulusSide(nanTrials)   = [];
        
        
        %% =================================================================
        % MODEL ACCURACY
        %% =================================================================
        
        modelAccuracy = 1 - double(modelStimulusSide == simulatedChoices);
        
        
        %% =================================================================
        % MODEL OUTPUT
        %
        % measureType = 1 --> RT proxy
        % measureType = 2 --> Confidence
        %% =================================================================
        
        if measureType == 1
            
            % RT proxy from confidence
            modelMeasure = 1 - confidenceAtDecision;
            
            % Align model mean with empirical mean
            modelMeasure = modelMeasure + ...
                (nanmean(behavioralMeasure) - ...
                 nanmean(nanmean(modelMeasure)));
            
        else
            
            % Confidence model output
            modelMeasure = modelOutput;
            
        end
        
        
        %% =================================================================
        % VECTORIZE MODEL VARIABLES
        %% =================================================================
        
        modelMeasure  = modelMeasure(:);
        modelAccuracy = modelAccuracy(:);
        modelEvidence = modelEvidence(:);
        
        
        %% =================================================================
        % EVIDENCE BINS (QUANTILES)
        %% =================================================================
        
        evidenceQuantiles = quantile(evidenceStrength(:,1), 3);
        
        [~, ~, evidenceBin_real] = histcounts( ...
            evidenceStrength(:,1), ...
            [-inf evidenceQuantiles inf]);
        
        [~, ~, evidenceBin_model] = histcounts( ...
            modelEvidence, ...
            [-inf evidenceQuantiles inf]);
        
        
        %% =================================================================
        % SPLIT CORRECT / ERROR TRIALS
        %% =================================================================
        
        correctTrials = find(trialAccuracy == 1);
        errorTrials   = find(trialAccuracy == 0);
        
        correctTrials_model = find(modelAccuracy == 1);
        errorTrials_model   = find(modelAccuracy == 0);
        
        
        %% #################################################################
        % CORRECT TRIALS
        %% #################################################################
        
        % Correlation between measure and evidence
        [pValue, correlationValue] = corr( ...
            zscore(behavioralMeasure(correctTrials,1)), ...
            zscore(evidenceStrength(correctTrials,1)), ...
            'type', 'spearman');
        
        correlations_correct(subject,:,measureType) = ...
            [pValue correlationValue];
        
        
        correctBins_real = evidenceBin_real(correctTrials);
        correctEvidence_real = evidenceStrength(correctTrials);
        correctMeasure_real  = behavioralMeasure(correctTrials);
        
        
        % Average per evidence bin
        for bin = 1:4
            
            meanMeasure_correct(subject,bin,measureType) = ...
                nanmean(correctMeasure_real(correctBins_real == bin,1));
            
            meanEvidence_correct(subject,bin,measureType) = ...
                nanmean(correctEvidence_real(correctBins_real == bin,1));
            
        end
        
        
        %% ---------------- MODEL : CORRECT TRIALS -------------------------
        
        correctBins_model = evidenceBin_model(correctTrials_model);
        correctEvidence_model = modelEvidence(correctTrials_model);
        correctMeasure_model  = modelMeasure(correctTrials_model);
        
        for bin = 1:4
            
            meanMeasure_correct_model(subject,bin,measureType) = ...
                nanmean(correctMeasure_model(correctBins_model == bin,1));
            
            meanEvidence_correct_model(subject,bin,measureType) = ...
                nanmean(correctEvidence_model(correctBins_model == bin,1));
            
        end
        
        
        %% #################################################################
        % ERROR TRIALS
        %% #################################################################
        
        [pValue, correlationValue] = corr( ...
            zscore(behavioralMeasure(errorTrials,1)), ...
            zscore(evidenceStrength(errorTrials,1)), ...
            'type', 'spearman');
        
        correlations_error(subject,:,measureType) = ...
            [pValue correlationValue];
        
        
        errorBins_real = evidenceBin_real(errorTrials);
        errorEvidence_real = evidenceStrength(errorTrials);
        errorMeasure_real  = behavioralMeasure(errorTrials);
        
        
        for bin = 1:4
            
            meanMeasure_error(subject,bin,measureType) = ...
                mean(errorMeasure_real(errorBins_real == bin,1));
            
            meanEvidence_error(subject,bin,measureType) = ...
                mean(errorEvidence_real(errorBins_real == bin,1));
            
        end
        
        
        %% ---------------- MODEL : ERROR TRIALS ---------------------------
        
        errorBins_model = evidenceBin_model(errorTrials_model);
        errorEvidence_model = modelEvidence(errorTrials_model);
        errorMeasure_model  = modelMeasure(errorTrials_model);
        
        for bin = 1:4
            
            meanMeasure_error_model(subject,bin,measureType) = ...
                mean(errorMeasure_model(errorBins_model == bin,1));
            
            meanEvidence_error_model(subject,bin,measureType) = ...
                mean(errorEvidence_model(errorBins_model == bin,1));
            
        end
        
    end
end


%% ========================================================================
% RESHAPE RESULTS FOR PLOTTING
%% ========================================================================

RT_correct_real  = reshape(meanMeasure_correct,        20, 4, 2);
EV_correct_real  = reshape(meanEvidence_correct,       20, 4, 2);

RT_error_real    = reshape(meanMeasure_error,          20, 4, 2);
EV_error_real    = reshape(meanEvidence_error,         20, 4, 2);

RT_correct_model = reshape(meanMeasure_correct_model,  20, 4, 2);
EV_correct_model = reshape(meanEvidence_correct_model, 20, 4, 2);

RT_error_model   = reshape(meanMeasure_error_model,    20, 4, 2);
EV_error_model   = reshape(meanEvidence_error_model,   20, 4, 2);

allCorr_correct = reshape(correlations_correct, 20, 2, 2);
allCorr_error   = reshape(correlations_error,   20, 2, 2);


%% ========================================================================
% FIGURE 1 - REACTION TIME
%% ========================================================================

figure(1)
clf


%% ------------------------------------------------------------------------
% PANEL 1 - RT vs Evidence
%% ------------------------------------------------------------------------

subplot(1,2,1)
hold on

% Correct trials - MODEL
shadedErrorBar( ...
    nanmean(EV_correct_model(:,:,1),1), ...
    nanmean(RT_correct_model(:,:,1),1), ...
    nanstd(RT_correct_model(:,:,1))./sqrt(20), ...
    'lineProps', {'-','color',[0.4660 0.6740 0.1880]});

% Correct trials - REAL
scatter( ...
    mean(EV_correct_real(:,:,1),1), ...
    mean(RT_correct_real(:,:,1),1), ...
    70, ...
    [0.4660 0.6740 0.1880], ...
    'filled');

% Error trials - MODEL
shadedErrorBar( ...
    nanmean(EV_error_model(:,:,1),1), ...
    nanmean(RT_error_model(:,:,1),1), ...
    nanstd(RT_error_model(:,:,1))./sqrt(20), ...
    'lineProps', {'-','color','r'});

% Error trials - REAL
scatter( ...
    mean(EV_error_real(:,:,1),1), ...
    mean(RT_error_real(:,:,1),1), ...
    70, ...
    'r', ...
    'filled', ...
    'MarkerEdgeColor',[1 1 1]);

xlabel('Evidence strength')
ylabel('Reaction time (s)')
set(gca,'fontsize',15)


%% ------------------------------------------------------------------------
% PANEL 2 - RT Correlation Coefficients
%% ------------------------------------------------------------------------

subplot(1,2,2)
hold on

% Group mean bars
bar(1, mean(allCorr_correct(:,1,1),1), ...
    'facecolor',[0.4660 0.6740 0.1880], ...
    'facealpha',0.5)

bar(2, mean(allCorr_error(:,1,1),1), ...
    'facecolor','r', ...
    'facealpha',0.5)

% Individual subjects
for subject = 1:20
    
    scatter( ...
        1 + (2*rand(1)-1)*0.05, ...
        allCorr_correct(subject,1,1), ...
        100, ...
        [0.4660 0.6740 0.1880], ...
        'filled', ...
        'markerfacealpha',0.5)
    
    scatter( ...
        2 + (2*rand(1)-1)*0.05, ...
        allCorr_error(subject,1,1), ...
        100, ...
        'r', ...
        'filled', ...
        'markerfacealpha',0.5)
    
end

xticks([])
yline(0, ...
    'linewidth',3, ...
    'color',[0.5 0.5 0.5], ...
    'linestyle','--')

ylabel('Scaling of RT with evidence')
set(gca,'fontsize',15)



%% ========================================================================
% FIGURE 2 - CONFIDENCE
%% ========================================================================

figure(2)
clf


%% ------------------------------------------------------------------------
% PANEL 1 - Confidence vs Evidence
%% ------------------------------------------------------------------------

subplot(1,2,1)
hold on

% Correct trials  MODEL
shadedErrorBar( ...
    nanmean(EV_correct_model(:,:,2),1), ...
    nanmean(RT_correct_model(:,:,2),1), ...
    nanstd(RT_correct_model(:,:,2))./sqrt(20), ...
    'lineProps', {'-','color',[0.4660 0.6740 0.1880]});

% Correct trials  REAL
scatter( ...
    mean(EV_correct_real(:,:,2),1), ...
    mean(RT_correct_real(:,:,2),1), ...
    70, ...
    [0.4660 0.6740 0.1880], ...
    'filled', ...
    'square', ...
    'MarkerEdgeColor',[1 1 1]);

% Error trials MODEL
shadedErrorBar( ...
    nanmean(EV_error_model(:,:,2),1), ...
    nanmean(RT_error_model(:,:,2),1), ...
    std(RT_error_model(:,:,2))./sqrt(20), ...
    'lineProps', {'-','color','r'});

% Error trials  REAL
scatter( ...
    mean(EV_error_real(:,:,2),1), ...
    mean(RT_error_real(:,:,2),1), ...
    70, ...
    'r', ...
    'filled', ...
    'MarkerEdgeColor',[1 1 1]);

xlabel('Evidence strength')
ylabel('Prob. High Confidence report')
set(gca,'fontsize',15)


%% ------------------------------------------------------------------------
% PANEL 2 - Confidence Correlation Coefficients
%% ------------------------------------------------------------------------

subplot(1,2,2)
hold on

% Group mean bars
bar(1, mean(allCorr_correct(:,1,2),1), ...
    'facecolor',[0.4660 0.6740 0.1880], ...
    'facealpha',0.5)

bar(2, mean(allCorr_error(:,1,2),1), ...
    'facecolor','r', ...
    'facealpha',0.5)

% Individual subjects
for subject = 1:20
    
    scatter( ...
        1 + (2*rand(1)-1)*0.05, ...
        allCorr_correct(subject,1,2), ...
        100, ...
        [0.4660 0.6740 0.1880], ...
        'filled', ...
        'markerfacealpha',0.5)
    
    scatter( ...
        2 + (2*rand(1)-1)*0.05, ...
        allCorr_error(subject,1,2), ...
        100, ...
        'r', ...
        'filled', ...
        'markerfacealpha',0.5)
    
end

xticks([])
yline(0, ...
    'linewidth',3, ...
    'color',[0.5 0.5 0.5], ...
    'linestyle','--')

ylabel('Scaling of Confidence report with evidence')
set(gca,'fontsize',15)