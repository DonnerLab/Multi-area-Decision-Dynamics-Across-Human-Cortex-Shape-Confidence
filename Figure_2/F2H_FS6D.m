%% =======================================================================
% Figure 2H and Figure S6D
%
% Across-trial correlation between the left and right DV at the end of
% stimulus viewing, after factoring out stimulus category and choice
% then averaged across the 4 categories
%
% 2H : neural DVs of PMd (left axis) and model DVs (right axis) of the
%      best-fitting model with the fitted feed-forward inhibition (FFI)
%      of each participant and with FFI = 0
% S6D: neural DVs of M1
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

% Single-hemisphere decoder output
decoderFolder = fullfile(projectFolder, 'data', 'decoder_output', 'PMd_M1_hemi');

% Path to functions
addpath(fullfile(projectFolder, 'functions'));

% Behavioral data of all subjects 
load(fullfile(behaviorFolder, 'data.mat'))


%% =======================================================================
% SETTINGS
%% =======================================================================

subjects = { ...
    'S01','S04','S06','S07','S12','S15','S17','S18','S11','S05', ...
    'S19','S20','S22','S10','S02','S08','S24','S21','S23','S25'};

nSubjects     = numel(subjects);
nPermutations = 10000;    % permutations of the permutation tests

% ROIs: name and areas whose decoder output is averaged
roiNames = {'PMd', 'M1'};
roiAreas = {{'6a','6d'}, {'4'}};

% Neural DVs: 0.9 - 1 s from stimulus onset (10 ms samples, first at -0.5 s)
timeSec       = ((1:200) - 1) * 0.01 - 0.5;
neuralWindow  = timeSec >= 0.9 & timeSec <= 1;

% Decoder of the model DVs: folds of the cross-validation within session
nFolds = 10;

%% =======================================================================
% NEURAL DVs (PMd: Figure 2H, M1: Figure S6D)
%% =======================================================================

% subjects x ROI
correlation_neural = nan(nSubjects,numel(roiNames));

for roi = 1:numel(roiNames)
    for subject = 1:nSubjects

        %% ---------------------------------------------------------------
        % LOAD THE DVs OF BOTH HEMISPHERES (all sessions)
        %% ---------------------------------------------------------------

        files = dir(fullfile(decoderFolder, ...
            [subjects{subject} '_SES*_l_CHOICEdecoder_' roiAreas{roi}{1} '_all.mat']));

        DV_right = [];   % left hemisphere
        DV_left  = [];   % right hemisphere
        choice   = [];   % -1 = left, +1 = right
        stimulus = [];   % -1 = left > right, +1 = right > left

        for file = 1:numel(files)

            session = sscanf(files(file).name, [subjects{subject} '_SES%d']);

            % DV of each hemisphere: decoder log odds in 0.9 - 1 s, averaged
            % over the areas of the ROI (trials x areas)
            areaDV_right = [];   % left hemisphere
            areaDV_left  = [];   % right hemisphere

            for area = 1:numel(roiAreas{roi})

                L = load(fullfile(decoderFolder, sprintf('%s_SES%d_l_CHOICEdecoder_%s_all.mat', ...
                    subjects{subject}, session, roiAreas{roi}{area})));
                R = load(fullfile(decoderFolder, sprintf('%s_SES%d_r_CHOICEdecoder_%s_all.mat', ...
                    subjects{subject}, session, roiAreas{roi}{area})));

                areaDV_right(:,area) = mean(double(L.log_odds_all(:,neuralWindow)),2);
                areaDV_left(:,area)  = mean(double(R.log_odds_all(:,neuralWindow)),2);

            end

            DV_right = [DV_right; mean(areaDV_right,2)];
            DV_left  = [DV_left;  mean(areaDV_left,2)];
            choice   = [choice;   double(L.Y_all(:))];
            stimulus = [stimulus; double(L.stimulus_all(:))];

        end


        %% ---------------------------------------------------------------
        % CORRELATION WITHIN EACH STIMULUS CATEGORY x CHOICE
        %% ---------------------------------------------------------------

        correlation_category = nan(2,2);

        for stimulusSide = [-1 1]
            for choiceSide = [-1 1]

                trials = stimulus == stimulusSide & choice == choiceSide;

                if sum(trials) >= 5
                    correlation_category((stimulusSide+3)/2,(choiceSide+3)/2) = ...
                        corr(DV_left(trials), DV_right(trials));
                end

            end
        end

        correlation_neural(subject,roi) = nanmean(correlation_category(:));

    end
end


%% =======================================================================
% MODEL DVs (Figure 2H)
%% =======================================================================

% subjects x FFI condition (1 = fitted FFI, 2 = FFI = 0)
correlation_model = nan(nSubjects,2);

for subject = 1:nSubjects

    %% -------------------------------------------------------------------
    % BEHAVIORAL DATA (the 3 drug conditions combined)
    %
    % columns: 1 = choice (0 = left, 1 = right), 2:11 = left contrasts,
    %          12:21 = right contrasts, 26 = MEG session
    %% -------------------------------------------------------------------

    subjectData = [];

    for drugCondition = 1:3
        subjectData = [subjectData; cell2mat(All_Data(subject,drugCondition))];
    end

    subjectData(isnan(subjectData(:,1)),:) = [];

    simulationInput.t     = 10;
    simulationInput.left  = subjectData(:,2:11)';
    simulationInput.right = subjectData(:,12:21)';

    nTrials           = size(subjectData,1);
    participantChoice = subjectData(:,1);
    stimulusSide      = double(mean(simulationInput.right - simulationInput.left, 1)' > 0);   % 1 = right > left
    sessionID         = subjectData(:,26);
    sessions          = unique(sessionID(~isnan(sessionID)))';

    load(fullfile(modelFolder, ['Parameters_model_9_sbj_' num2str(subject) '.mat']));
    bestParameters = Parameters(subject,:,2);

    for ffiCondition = 1:2

        parameters = bestParameters;

        if ffiCondition == 2
            parameters(5) = 0;   % no feed-forward inhibition
        end

        %% ---------------------------------------------------------------
        % SIMULATION (one per trial) WITH ABSORBING BOUND
        %
        % P_left, P_right: 100 time steps of 10 ms (stimulus) x trials;
        % winning DV stays at the bound from its crossing on, losing DV
        % accumulates until it reaches the bound
        %% ---------------------------------------------------------------

        [~, Choices, ~, ~, bound, ~, ~, P_left, P_right, ~, ~, Cross_step] = ...
            leaky_integrator_simulations_fine(parameters, simulationInput, 4, 1);

        nSteps = size(P_left,1);

        % apply absorbing bound 
        [P_left, P_right] = absorb_bound_dvs( ...
            reshape(P_left,nSteps,nTrials), reshape(P_right,nSteps,nTrials), ...
            bound, Choices(:) == 1, Cross_step(:));

        % Model DVs at the end of stimulus viewing (last time step, 1 s)
        DV_right = P_right(end,:)';
        DV_left  = P_left(end,:)';

        % Model choice (decoder labels): +1 = right, -1 = left
        modelChoice = double(Choices(:));
        modelChoice(modelChoice == 0) = -1;


        %% ---------------------------------------------------------------
        % DECODER LOG ODDS, as for the neural DVs
        %
        % Within each session, 10-fold cross-validation: right DV ->
        % right choice, left DV -> left choice; classes balanced by
        % random downsampling of the larger class; DV z-scored with the
        % training mean and SD; logistic regression (fitclinear)
        %% ---------------------------------------------------------------

        logOdds_right = nan(nTrials,1);
        logOdds_left  = nan(nTrials,1);

        for session = sessions

            sessionTrials = find(sessionID == session);

            if numel(sessionTrials) < nFolds || numel(unique(modelChoice(sessionTrials))) < 2
                continue
            end

            folds = cvpartition(numel(sessionTrials), 'KFold', nFolds);

            for side = 1:2

                if side == 1
                    X_all = DV_right(sessionTrials);
                    Y_all = modelChoice(sessionTrials);
                else
                    X_all = DV_left(sessionTrials);
                    Y_all = -modelChoice(sessionTrials);
                end

                sessionLogOdds = nan(numel(sessionTrials),1);

                for fold = 1:nFolds

                    X_train = X_all(training(folds,fold));
                    Y_train = Y_all(training(folds,fold));

                    % Balance the classes (downsample the larger class)
                    negative = find(Y_train == -1);
                    positive = find(Y_train ==  1);

                    if isempty(negative) || isempty(positive)
                        continue
                    end

                    if numel(negative) > numel(positive)
                        balanced = [positive; negative(randperm(numel(negative), numel(positive)))];
                    else
                        balanced = [negative; positive(randperm(numel(positive), numel(negative)))];
                    end

                    X_train = X_train(balanced);
                    Y_train = Y_train(balanced);

                    % z-score with the training mean and SD, logistic regression
                    mu    = mean(X_train);
                    sigma = std(X_train);
                    if sigma == 0, sigma = 1; end

                    decoder = fitclinear((X_train - mu) ./ sigma, Y_train, 'Learner', 'logistic');

                    % Log odds of the test trials
                    testTrials = test(folds,fold);
                    sessionLogOdds(testTrials) = (X_all(testTrials) - mu) ./ sigma .* decoder.Beta + decoder.Bias;

                end

                if side == 1
                    logOdds_right(sessionTrials) = sessionLogOdds;
                else
                    logOdds_left(sessionTrials)  = sessionLogOdds;
                end

            end

        end


        %% ---------------------------------------------------------------
        % CORRELATION WITHIN EACH STIMULUS CATEGORY x CHOICE 
        %% ---------------------------------------------------------------

        correlation_category = nan(2,2);

        for stimulus = 0:1
            for choice = 0:1

                trials = stimulusSide == stimulus & participantChoice == choice & ...
                    isfinite(logOdds_right) & isfinite(logOdds_left);

                if sum(trials) >= 5
                    correlation_category(stimulus+1,choice+1) = ...
                        corr(logOdds_right(trials), logOdds_left(trials));
                end

            end
        end

        % Mean of the 4 categories
        correlation_model(subject,ffiCondition) = nanmean(correlation_category(:));

    end

end


%% =======================================================================
% STATISTICS

%% =======================================================================

pValue_PMd = permutationTest(correlation_neural(:,1), zeros(nSubjects,1), nPermutations);
pValue_M1  = permutationTest(correlation_neural(:,2), zeros(nSubjects,1), nPermutations);

pValue_fittedFFI = permutationTest(correlation_model(:,1), zeros(nSubjects,1), nPermutations);
pValue_zeroFFI   = permutationTest(correlation_model(:,2), zeros(nSubjects,1), nPermutations);
pValue_FFI_diff  = permutationTest(correlation_model(:,1) - correlation_model(:,2), ...
    zeros(nSubjects,1), nPermutations);


%% =======================================================================
% FIGURE 1 (Figure 2H) - PMd NEURAL DVs AND MODEL DVs
%
% Neural DVs on the left y-axis, model DVs on the right y-axis; the right
% axis is the left axis scaled by modelScale (zeros aligned)
%% =======================================================================

modelScale = max(abs(correlation_model(:))) / max(abs(correlation_neural(:,1)));

figure(1)

yyaxis left
hold on

bar(1, mean(correlation_neural(:,1)), 0.6, 'FaceColor', [0.80 0.80 0.80], 'EdgeColor', 'k');
scatter(1 + (rand(nSubjects,1)-0.5)*0.25, correlation_neural(:,1), 50, 'k', 'filled', 'MarkerFaceAlpha', 0.5);

bar(2, mean(correlation_model(:,1)) / modelScale, 0.6, 'FaceColor', [0.96 0.67 0.55], 'EdgeColor', 'k');
scatter(2 + (rand(nSubjects,1)-0.5)*0.25, correlation_model(:,1) / modelScale, 50, 'k', 'filled', 'MarkerFaceAlpha', 0.5);

bar(3, mean(correlation_model(:,2)) / modelScale, 0.6, 'FaceColor', [0.99 0.85 0.70], 'EdgeColor', 'k');
scatter(3 + (rand(nSubjects,1)-0.5)*0.25, correlation_model(:,2) / modelScale, 50, 'k', 'filled', 'MarkerFaceAlpha', 0.5);

plot([0.4 3.6], [0 0], 'k-');

leftLimits = ylim;
ylabel({'Correlation coefficient', '(neural DVs)'})
set(gca,'ycolor','k')

yyaxis right
ylim(leftLimits * modelScale)
ylabel({'Correlation coefficient', '(model DVs)'})
set(gca,'ycolor','k')

xlim([0.4 3.6])
xticks(1:3)
xticklabels({'Neural DVs', 'Model DVs', 'Model DVs (FFI = 0)'})
title({sprintf('neural: p = %.2g,  model: p = %.2g,  FFI = 0: p = %.2g', ...
    pValue_PMd, pValue_fittedFFI, pValue_zeroFFI); ...
    sprintf('model vs FFI = 0: p = %.2g', pValue_FFI_diff)}, 'FontWeight', 'normal')
box off
set(gca,'FontSize',12)


%% =======================================================================
% FIGURE 2 (Figure S6D) - M1 NEURAL DVs
%% =======================================================================

figure(2)
hold on

bar(1, mean(correlation_neural(:,2)), 'FaceColor', 'k', 'FaceAlpha', 0.5)

for subject = 1:nSubjects
    scatter(1 + (rand(1)*2 - 1)*0.05, correlation_neural(subject,2), 100, 'k', ...
        'filled', 'MarkerFaceAlpha', 0.5)
end

yline(0,'k-');
xticks([])
ylabel({'Correlation coefficient', '(neural DVs)'})
title(sprintf('M1: p = %.2g', pValue_M1), 'FontWeight', 'normal')
box off
set(gca,'FontSize',15)


