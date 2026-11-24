%% =======================================================================
% Figure S9D-E
%
% Regression of the model decision time on evidence strength and on the
% winning and losing model DV, over time, for the best-fitting model
% (model 9 parameters) with
%
%   collapsing bound (simulation 4, solid lines)
%   fixed bound      (simulation 2, dashed lines)
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
nSimulations  = 1000;       % simulations per trial
nPermutations = 10000;    % permutations of the permutation tests

% Bounds: 4 = collapsing bound, 2 = fixed bound
boundNames       = {'Collapsing bound', 'Fixed bound'};
simulationNumber = [4, 2];

% Model time steps: 100 steps of 10 ms during the stimulus
nSteps = 100;

% Fixed bound: every trial must cross within 500 ms (= step 50)
crossingLimit = 50;
nSimulations_calibration = 5;
nIterations_calibration  = 18;

% Non-overlapping windows of 5 steps (50 ms); window centers (s)
windowLength = 5;
windowStarts = 1:windowLength:(nSteps - windowLength + 1);
nWindows     = numel(windowStarts);
windowTime   = ((windowStarts - 1 + windowLength/2) * 10) / 1000;

% End of stimulus viewing: windows in the last 100 ms
offsetWindows = find(windowTime > 0.9);

% Decoder of the model DVs: folds of the cross-validation within session
nFolds = 10;

% Ridge penalty of the model regression
lambda = 1e-3;


% subjects x windows x [winning DV, losing DV] x bound
weights = nan(nSubjects,nWindows,2,2);


%% =======================================================================
% SIMULATIONS AND REGRESSIONS
%% =======================================================================

for bound = 1:2
    for subject = 1:nSubjects

        %% ---------------------------------------------------------------
        % STIMULI AND PARAMETERS
        %% ---------------------------------------------------------------

        subjectData = [];

        for drugCondition = 1:3
            subjectData = [subjectData; cell2mat(All_Data(subject,drugCondition))];
        end

        subjectData(isnan(subjectData(:,1)),:) = [];

        simulationInput.t     = 10;
        simulationInput.left  = subjectData(:,2:11)';
        simulationInput.right = subjectData(:,12:21)';

        nTrials   = size(subjectData,1);
        sessionID = subjectData(:,26);

        load(fullfile(modelFolder, ['Parameters_model_9_sbj_' num2str(subject) '.mat']));
        parameters = Parameters(subject,:,2);


        %% ---------------------------------------------------------------
        % FIXED BOUND: HIGHEST LEVEL AT WHICH ALL TRIALS CROSS IN 500 ms
        %
        % Bisection on a log scale between 1e-6 x and 1 x the initial
        % bound of the model (parameter 8, x 15)
        %% ---------------------------------------------------------------

        if simulationNumber(bound) == 2

            initialBound = abs(parameters(8)) * 15;
            lowerLevel   = log(1e-6 * initialBound);
            upperLevel   = log(initialBound);

            for iteration = 1:nIterations_calibration

                testLevel = 0.5 * (lowerLevel + upperLevel);
                testParameters = parameters;
                testParameters(8) = exp(testLevel) / 15;
                testParameters(11) =0;

                [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, Cross_step] = leaky_integrator_simulations_fine( ...
                    testParameters, simulationInput, 2, nSimulations_calibration);

                if all(Cross_step(:) <= crossingLimit)
                    lowerLevel = testLevel;
                else
                    upperLevel = testLevel;
                end

            end

            parameters(8) = exp(lowerLevel) / 15;
            parameters(11)=0;
        end


        %% ---------------------------------------------------------------
        % SIMULATION (nSimulations per trial)
        %
        % Fixed bound: lowered by 10% until all simulated trials cross
        % within 500 ms
        %% ---------------------------------------------------------------

        [~, Choices, ~, ~, boundValues, ~, ~, P_left, P_right, ~, Cross_time, Cross_step] = ...
            leaky_integrator_simulations_fine(parameters, simulationInput, simulationNumber(bound), nSimulations);

        nShrink = 0;

        while simulationNumber(bound) == 2 && ~all(Cross_step(:) <= crossingLimit) && nShrink < 60

            parameters(8) = 0.9 * parameters(8);
            parameters(11) =0;
            nShrink = nShrink + 1;

            [~, Choices, ~, ~, boundValues, ~, ~, P_left, P_right, ~, Cross_time, Cross_step] = ...
                leaky_integrator_simulations_fine(parameters, simulationInput, simulationNumber(bound), nSimulations);

        end

        % One row per simulated trial (trials x simulations)
        isRight      = Choices(:) == 1;
        crossed      = ~isnan(Cross_step(:));
        decisionTime = Cross_time(:);
        evidence     = repmat(abs(mean(simulationInput.right - simulationInput.left,1))', nSimulations, 1);
        session      = repmat(sessionID, nSimulations, 1);
        sessions     = unique(session(~isnan(session)))';

        % Model choice (decoder labels): +1 = right, -1 = left
        modelChoice = double(isRight);
        modelChoice(~isRight) = -1;

        % Accumulators: simulated trials x time steps, absorbing bound
        accumulator_left  = reshape(permute(P_left,  [2 3 1]), nTrials*nSimulations, nSteps);
        accumulator_right = reshape(permute(P_right, [2 3 1]), nTrials*nSimulations, nSteps);

        [accumulator_left, accumulator_right] = absorb_bound_dvs( ...
            accumulator_left', accumulator_right', boundValues, isRight, Cross_step(:));

        accumulator_left  = accumulator_left';
        accumulator_right = accumulator_right';


        %% ---------------------------------------------------------------
        % IN EACH 50 ms WINDOW: DECODER LOG ODDS AND REGRESSION
        %% ---------------------------------------------------------------

        for window = 1:nWindows

            windowSteps = windowStarts(window) + (0:windowLength-1);

            DV_left  = mean(accumulator_left(:,windowSteps),2);
            DV_right = mean(accumulator_right(:,windowSteps),2);

            % Decoder log odds, as for the neural DVs: within each session,
            % 10-fold cross-validation; right DV -> right choice, left DV
            % -> left choice; smaller class upsampled in the training set;
            % DV z-scored with the training mean and SD; logistic
            % regression (fitclinear)
            logOdds_right = nan(nTrials*nSimulations,1);
            logOdds_left  = nan(nTrials*nSimulations,1);

            for thisSession = sessions

                sessionTrials = find(session == thisSession);
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

                        % Upsample the smaller class
                        negative = find(Y_train == -1);
                        positive = find(Y_train ==  1);

                        if isempty(negative) || isempty(positive)
                            continue
                        end

                        if numel(negative) > numel(positive)
                            extra = positive(randi(numel(positive), numel(negative) - numel(positive), 1));
                        else
                            extra = negative(randi(numel(negative), numel(positive) - numel(negative), 1));
                        end

                        X_train = [X_train; X_train(extra)];
                        Y_train = [Y_train; Y_train(extra)];

                        % z-score with the training mean and SD, logistic regression
                        [X_train, mu, sigma] = zscore(X_train);
                        decoder = fitclinear(X_train, Y_train, 'Learner', 'logistic');

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

            % Regression (trials with a bound crossing):
            % decision time ~ evidence + winning DV + losing DV
            winningDV = logOdds_left;   winningDV(isRight) = logOdds_right(isRight);
            losingDV  = logOdds_right;  losingDV(isRight)  = logOdds_left(isRight);

            trials = crossed & isfinite(winningDV) & isfinite(losingDV);

            if sum(trials)>0
            regressors   = [evidence(trials), winningDV(trials), losingDV(trials)];
            trialTime    = decisionTime(trials);
            trialSession = session(trials);

            % z-scored within each session
            for thisSession = unique(trialSession)'
                inSession = trialSession == thisSession;
                regressors(inSession,:) = zscore(regressors(inSession,:));
                trialTime(inSession)    = zscore(trialTime(inSession));
            end

            
            rawWinningDV = DV_left;   rawWinningDV(isRight) = DV_right(isRight);
            rawLosingDV  = DV_right;  rawLosingDV(isRight)  = DV_left(isRight);



            RTmodel = fitrlinear(regressors, trialTime, ...
                'Learner', 'leastsquares', 'Regularization', 'ridge', 'Lambda', lambda);

            % Weights of the winning and losing DV
            weights(subject,window,:,bound) = RTmodel.Beta(2:3);
            else
           weights(subject,window,:,bound) = [];
    
            end
        end

    end
end


%% =======================================================================
% |WINNING| - |LOSING| WEIGHT AT THE END OF STIMULUS VIEWING
%% =======================================================================

% subjects x bound
weightDifference = nan(nSubjects,2);

for bound = 1:2
    offsetWeights = squeeze(nanmean(weights(:,offsetWindows,:,bound),2));
    weightDifference(:,bound) = abs(offsetWeights(:,1)) - abs(offsetWeights(:,2));
end

pValues = [permutationTest(weightDifference(:,1), zeros(nSubjects,1), nPermutations), ...
           permutationTest(weightDifference(:,2), zeros(nSubjects,1), nPermutations)];


%% =======================================================================
% FIGURE 1 (Figure S9D) - WEIGHTS OVER TIME
%
% black: winning DV, blue: losing DV; solid: collapsing, dashed: fixed
%% =======================================================================

dvColors   = [0 0 0; 0.15 0.45 0.85];
lineStyles = {'-', '--'};

figure(1)
hold on

for bound = 1:2
    for dv = 1:2

        shadedErrorBar( ...
            windowTime, ...
            nanmean(weights(:,:,dv,bound),1), ...
            nanstd(weights(:,:,dv,bound),0,1) ./ sqrt(nSubjects), ...
            'lineProps', {lineStyles{bound},'color',dvColors(dv,:),'LineWidth',2.5});

    end
end

yline(0, 'Color', [0.3 0.3 0.3]);
xlim([0 1])
xlabel('Time from stimulus onset (s)')
ylabel('Weight (a.u.)')
title('Winning (black) and losing (blue) DV; collapsing (solid), fixed bound (dashed)', ...
    'FontWeight', 'normal', 'FontSize', 10)
box off
set(gca, 'FontSize', 12)


%% =======================================================================
% FIGURE 2 (Figure S9E) - |WINNING| - |LOSING| WEIGHT AT STIMULUS OFFSET
%
% Fixed bound on the right y-axis: left y-axis scaled by boundScale
% (zeros aligned)
%% =======================================================================

boundScale = max(abs(weightDifference(:,2))) / max(abs(weightDifference(:,1)));

figure(2)

yyaxis left
hold on

bar(1, nanmean(weightDifference(:,1)), 0.7, 'FaceColor', [0.96 0.67 0.55], 'EdgeColor', 'k');
scatter(1 + (rand(nSubjects,1) - 0.5)*0.3, weightDifference(:,1), 70, 'k', 'filled', 'MarkerFaceAlpha', 0.5);

bar(2, nanmean(weightDifference(:,2)) / boundScale, 0.7, 'FaceColor', [0.96 0.67 0.55], 'EdgeColor', 'k');
scatter(2 + (rand(nSubjects,1) - 0.5)*0.3, weightDifference(:,2) / boundScale, 70, 'k', 'filled', 'MarkerFaceAlpha', 0.5);

plot([0.3 2.7], [0 0], 'k-');

leftLimits = ylim;
ylabel({'Difference in absolute values', 'of winning and losing DV weights'})
set(gca,'ycolor','k')

yyaxis right
ylim(leftLimits * boundScale)
set(gca,'ycolor','k')

xlim([0.3 2.7])
xticks(1:2)
xticklabels(boundNames)
title(sprintf('collapsing: p = %.2g,  fixed: p = %.2g', pValues(1), pValues(2)), 'FontWeight', 'normal')
box off
set(gca, 'FontSize', 12)
