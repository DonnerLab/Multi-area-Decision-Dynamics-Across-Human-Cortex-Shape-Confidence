%% =======================================================================
% Figure S6E
%
% Within-trial correlation between the winning and losing model DV in the
% early stimulus presentation (0 - 0.5 s), shown as cross-correlograms,
% for the best-fitting model with feed-forward inhibition (FFI) = 0 and
% with the fitted FFI of each participant
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

nSubjects = 20;

% Early stimulus presentation: first 51 time steps of the model (10 ms each)
earlyWindow = 1:51;

% Cross-correlation lags: +-30 time steps = +-300 ms
maxLag  = 30;
lagAxis = (-maxLag:maxLag) * 10;

% subjects x lags x FFI condition (1 = FFI = 0, 2 = fitted FFI)
crossCorrelation = nan(nSubjects, 2*maxLag+1, 2);



%% =======================================================================
% MAIN LOOP ACROSS SUBJECTS
%% =======================================================================

for subject = 1:nSubjects

    % Stimulus of each trial (the 3 drug conditions combined)
    subjectData = [];

    for drugCondition = 1:3
        subjectData = [subjectData; cell2mat(All_Data(subject,drugCondition))];
    end

    subjectData(isnan(subjectData(:,1)),:) = [];

    simulationInput.t     = 10;
    simulationInput.left  = subjectData(:,2:11)';
    simulationInput.right = subjectData(:,12:21)';

    load(fullfile(modelFolder, ['Parameters_model_9_sbj_' num2str(subject) '.mat']));
    bestParameters = Parameters(subject,:,2);

    for ffiCondition = 1:2

        parameters = bestParameters;

        if ffiCondition == 1
            parameters(5) = 0;   % no feed-forward inhibition
        end

        %% ---------------------------------------------------------------
        % SIMULATION (one per trial): accumulators, trials x time
        %% ---------------------------------------------------------------

        [~, Choices, ~, ~, ~, ~, ~, P_left, P_right] = ...
            leaky_integrator_simulations_fine(parameters, simulationInput, 4, 1);

        nTimeSteps  = size(P_left,1);
        DV_left     = reshape(P_left,  nTimeSteps, [])';
        DV_right    = reshape(P_right, nTimeSteps, [])';
        rightChoice = Choices(:) == 1;
        nTrials     = numel(rightChoice);


        %% ---------------------------------------------------------------
        % DECODER LOG ODDS AT EACH TIME POINT AND WINNING / LOSING DV
        %% ---------------------------------------------------------------

        logOdds_right = nan(nTrials,nTimeSteps);
        logOdds_left  = nan(nTrials,nTimeSteps);

        for time = earlyWindow
            [logOdds_right(:,time), logOdds_left(:,time)] = model_logodds( ...
                DV_right, DV_left, rightChoice, ones(nTrials,1), time);
        end

        winningDV = logOdds_left;
        losingDV  = logOdds_right;
        winningDV(rightChoice,:) = logOdds_right(rightChoice,:);
        losingDV(rightChoice,:)  = logOdds_left(rightChoice,:);


        %% ---------------------------------------------------------------
        % SINGLE-TRIAL CROSS-CORRELATION (0 - 0.5 s)
        %
 
        %% ---------------------------------------------------------------

        trialCorrelation = nan(nTrials, 2*maxLag+1);

        for trial = 1:nTrials

            winning = winningDV(trial,earlyWindow);
            losing  = losingDV(trial,earlyWindow);

            valid = ~isnan(winning) & ~isnan(losing);

            if sum(valid) < 3
                continue
            end

            winning = zscore(detrend(winning(valid),2));
            losing  = zscore(detrend(losing(valid),2));

            if ~any(isnan(winning)) && ~any(isnan(losing))
                trialCorrelation(trial,:) = xcorr(winning', losing', maxLag, 'normalized');
            end

        end

        crossCorrelation(subject,:,ffiCondition) = nanmean(trialCorrelation,1);

    end

end


%% =======================================================================
% FIGURE (Figure S6E)
%% =======================================================================

ffiColors = [1 0 0; 1 1 0];   % red: FFI = 0, yellow: fitted FFI

figure
hold on

for ffiCondition = 1:2
    shadedErrorBar( ...
        lagAxis, ...
        nanmean(crossCorrelation(:,:,ffiCondition),1), ...
        nanstd(crossCorrelation(:,:,ffiCondition),0,1) ./ sqrt(nSubjects), ...
        'lineProps', {'-','color',ffiColors(ffiCondition,:),'LineWidth',1.5});
end

xline(0,'k--');
yline(0,'k:');
xlabel('Time lag (ms)')
ylabel('Correlation coefficient')
title('Inhibition: 0 (red), fitted (yellow)', 'FontWeight', 'normal')
box off
set(gca,'FontSize',15)
