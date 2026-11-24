%% =======================================================================
% Figure 3D and Figure S8B
%
% Confidence ~ evidence + winning DV + losing DV at the end of stimulus
% viewing (logistic regression): |weight of winning DV| - |weight of
% losing DV| for neural and model DVs
%
% 3D : PMd neural DVs and DVs of the best-fitting model (model 9,
%      collapsing bound, confidence asymmetry weight)
% S8B: PMd neural DVs and DVs of the model with rectification of the DVs
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

% Single-hemisphere decoder output with confidence and RT of each trial:
% <subject>_SES<session>_<l/r>_CHOICEdecoder_<area>_all.mat
decoderFolder = fullfile(projectFolder, 'data', 'decoder_output', 'PMd_M1_hemi_behavior');

% Path to functions
addpath(fullfile(projectFolder, 'functions'));

% Behavioral data of all subjects (All_Data: subjects x drug conditions)
load(fullfile(behaviorFolder, 'data.mat'))


%% =======================================================================
% SETTINGS
%% =======================================================================

subjects = { ...
    'S01','S02','S04','S05','S06','S07','S08','S10','S11','S12', ...
    'S15','S17','S18','S19','S20','S21','S22','S23','S24','S25'};

nSubjects     = numel(subjects);
nPermutations = 10000;    % permutations of the permutation tests

% End of stimulus viewing (samples of the decoder time base)
offsetWindow = 134:168;

% Model DVs: 100 time steps of 10 ms during the stimulus (step k = DV
% after k x 10 ms); end of stimulus viewing: 0.9 - 1 s
modelWindow = 90:100;

% Decoder of the model DVs: folds of the cross-validation within session
nFolds = 10;

% Ridge penalty of the model regression
lambda = 1e-3;


% Models: Figure 3D, Figure S8B
modelNames       = {'Best-fitting model (Figure 3D)', 'Model with rectification (Figure S8B)'};
modelNumbers     = [9, 14];
simulationNumber = [4, 5];   % 4 = collapsing bound, 5 = collapsing bound + rectification


%% =======================================================================
% NEURAL DVs (PMd)
%% =======================================================================

betas_neural = nan(nSubjects,3);

for subject = 1:nSubjects

    files = dir(fullfile(decoderFolder, [subjects{subject} '_SES*_l_CHOICEdecoder_6A6D_all.mat']));

    DV_R       = [];   % left hemisphere, DV for the right choice
    DV_L       = [];   % right hemisphere, DV for the left choice
    isRight    = [];
    confidence = [];
    RT         = [];
    evidence   = [];
    session    = [];

    for file = 1:numel(files)

        sessionNumber = sscanf(files(file).name, [subjects{subject} '_SES%d']);

        L = load(fullfile(decoderFolder, sprintf('%s_SES%d_l_CHOICEdecoder_6A6D_all.mat', subjects{subject}, sessionNumber)));
        R = load(fullfile(decoderFolder, sprintf('%s_SES%d_r_CHOICEdecoder_6A6D_all.mat', subjects{subject}, sessionNumber)));

        DV_R       = [DV_R;       double(L.log_odds_all(:,1:201))];
        DV_L       = [DV_L;       double(R.log_odds_all(:,1:201))];
        isRight    = [isRight;    double(L.Y_all(:)) == 1];
        confidence = [confidence; double(L.confidence_all(:)) - 1];
        RT         = [RT;         double(L.rts_all(:))];
        evidence   = [evidence;   double(L.evidence_strength(:))];
        session    = [session;    sessionNumber * ones(size(L.log_odds_all,1),1)];

    end

    % Trials with a response
    keep    = RT > 0;
    isRight = logical(isRight);

    % Winning and losing DV (mean over the samples of the window)
    DV_right = mean(DV_R(keep,offsetWindow),2);
    DV_left  = mean(DV_L(keep,offsetWindow),2);

    winningDV = DV_left;   winningDV(isRight(keep)) = DV_right(isRight(keep));
    losingDV  = DV_right;  losingDV(isRight(keep))  = DV_left(isRight(keep));

    evidence = evidence(keep);
    session  = session(keep);

    % z-scored within each session
    for thisSession = unique(session)'
        trials = session == thisSession;
        evidence(trials)  = zscore(evidence(trials));
        winningDV(trials) = zscore(winningDV(trials));
        losingDV(trials)  = zscore(losingDV(trials));
    end

    % confidence ~ evidence + winning DV + losing DV
    confidenceModel = fitglm([evidence, winningDV, losingDV], confidence(keep), 'Distribution', 'binomial');

    betas_neural(subject,:) = confidenceModel.Coefficients.Estimate(2:4)';

end

asymmetry_neural = abs(betas_neural(:,2)) - abs(betas_neural(:,3));


%% =======================================================================
% MODEL DVs
%% =======================================================================

% models x subjects
asymmetry_model = nan(numel(modelNumbers),nSubjects);

for model = 1:numel(modelNumbers)
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

        nTrials   = size(subjectData,1);
        sessionID = subjectData(:,26);
        sessions  = unique(sessionID(~isnan(sessionID)))';

        % Fitted parameters; model with rectification: inhibition = 1
        load(fullfile(modelFolder, ['Parameters_model_' num2str(modelNumbers(model)) ...
            '_sbj_' num2str(subject) '.mat']));
        parameters = Parameters(subject,:,2);

        if modelNumbers(model) == 14
            parameters(5) = 1;
        end


        %% ---------------------------------------------------------------
        % SIMULATION (one per trial) WITH ABSORBING BOUND
        %
        %% ---------------------------------------------------------------

        [~, Choices, ~, modelConfidence, bound, ~, ~, P_left, P_right, ~, ~, Cross_step] = ...
            leaky_integrator_simulations_fine(parameters, simulationInput, simulationNumber(model), 1);

        nSteps = size(P_left,1);

        [P_left, P_right] = absorb_bound_dvs( ...
            reshape(P_left,nSteps,nTrials), reshape(P_right,nSteps,nTrials), ...
            bound, Choices(:) == 1, Cross_step(:));

        modelRight      = Choices(:) == 1;
        modelConfidence = double(modelConfidence(:));     % 1 = high, 0 = low
        modelEvidence   = abs(mean(simulationInput.right - simulationInput.left,1))';

        % Model choice (decoder labels): +1 = right, -1 = left
        modelChoice = double(modelRight);
        modelChoice(modelChoice == 0) = -1;

        % Model DVs at the end of stimulus viewing (0.9 - 1 s)
        DV_right = mean(P_right(modelWindow,:),1)';
        DV_left  = mean(P_left(modelWindow,:),1)';


        %% ---------------------------------------------------------------
        % DECODER LOG ODDS, as for the neural DVs
        %
  
        %% ---------------------------------------------------------------

        logOdds_right = nan(nTrials,1);
        logOdds_left  = nan(nTrials,1);

        for session = sessions

            sessionTrials = find(sessionID == session);
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


        %% ---------------------------------------------------------------
        % REGRESSION: model confidence ~ evidence + winning DV + losing DV
        %% ---------------------------------------------------------------

        valid = isfinite(logOdds_right) & isfinite(logOdds_left);

        winningDV = logOdds_left;   winningDV(modelRight) = logOdds_right(modelRight);
        losingDV  = logOdds_right;  losingDV(modelRight)  = logOdds_left(modelRight);

        winningDV    = winningDV(valid);
        losingDV     = losingDV(valid);
        evidence     = modelEvidence(valid);
        trialSession = sessionID(valid);

        % z-scored within each session
        for thisSession = unique(trialSession)'
            trials = trialSession == thisSession;
            evidence(trials)  = zscore(evidence(trials));
            winningDV(trials) = zscore(winningDV(trials));
            losingDV(trials)  = zscore(losingDV(trials));
        end

        confidenceModel = fitclinear([evidence, winningDV, losingDV], modelConfidence(valid), ...
            'Learner', 'logistic', 'Regularization', 'ridge', 'Lambda', lambda);

        betas_model = confidenceModel.Beta';

        asymmetry_model(model,subject) = abs(betas_model(2)) - abs(betas_model(3));

    end
end


%% =======================================================================
% FIGURE (Figure 3D, Figure S8B)
%

%% =======================================================================

pValue_neural = permutationTest(asymmetry_neural, zeros(nSubjects,1), nPermutations);

figure('Color','w','Position',[200 200 700 380]);

for model = 1:numel(modelNumbers)

    pValue_model = permutationTest(asymmetry_model(model,:)', zeros(nSubjects,1), nPermutations);

    modelScale = max(abs(asymmetry_model(model,:))) / max(abs(asymmetry_neural));

    subplot(1,numel(modelNumbers),model)

    yyaxis left
    hold on

    bar(1, nanmean(asymmetry_neural), 0.6, 'FaceColor', [0.80 0.80 0.80], 'EdgeColor', 'k');
    scatter(1 + (rand(nSubjects,1)-0.5)*0.25, asymmetry_neural, 50, 'k', 'filled', 'MarkerFaceAlpha', 0.5);

    bar(2, nanmean(asymmetry_model(model,:)) / modelScale, 0.6, 'FaceColor', [0.96 0.67 0.55], 'EdgeColor', 'k');
    scatter(2 + (rand(nSubjects,1)-0.5)*0.25, asymmetry_model(model,:) / modelScale, 50, 'k', 'filled', 'MarkerFaceAlpha', 0.5);

    plot([0.4 2.6], [0 0], 'k-');

    leftLimits = ylim;
    ylabel({'Difference in absolute values', 'of winning and losing DV weights'})
    set(gca,'ycolor','k')

    yyaxis right
    ylim(leftLimits * modelScale)
    set(gca,'ycolor','k')

    xlim([0.4 2.6])
    xticks([1 2])
    xticklabels({'Neural DVs','Model DVs'})
    title({modelNames{model}; sprintf('neural: p = %.2g,  model: p = %.2g', pValue_neural, pValue_model)}, ...
        'FontWeight', 'normal', 'FontSize', 10)
    box off
    set(gca,'FontSize',12)


   
end