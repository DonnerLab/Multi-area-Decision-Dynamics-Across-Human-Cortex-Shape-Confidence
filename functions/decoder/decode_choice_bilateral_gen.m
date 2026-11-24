function decode_choice_bilateral_gen(subject, sessions, Rec)
%% =======================================================================
% Temporal generalization (train time x test time) of the bilateral
% choice decoder (Figure 2C: PMd, Figure S4B: M1)
%
% One decoder per area group (vertices of both hemispheres pooled),
% trained on all sessions together:
%   PMd = areas 6a + 6d
%   M1  = area 4
%

% Inputs:
%   subject  - subject number (e.g. 1 for S01)
%   sessions - session number(s)
%   Rec      - recording number(s) of the session
%
% Output (one file per area group):
%   data/generalization_matrix/
%     <subject>_<sessions>_<PMd|M1>_CHOICEdecoder_GENERALIZATION.mat
%     accuracy_gen: training time (rows) x testing time (columns),
%     mean over folds
%
%% =======================================================================


%% =======================================================================
% PATHS
%% =======================================================================

erfFolder    = '~/Desktop/2AFC/ERF';
trialsFolder = '~/Desktop/2AFC/trials';

functionFolder = fileparts(mfilename('fullpath'));
outputFolder   = fullfile(fileparts(fileparts(functionFolder)), 'data', 'generalization_matrix');


%% =======================================================================
% SETTINGS
%% =======================================================================

rng(42);

% Subject ID; participants with an additional session have session 1
% stored as session 11
subjectNumber = subject;
subject = sprintf('S%02d', subjectNumber);

if sessions == 1 && any([1,2,4,5,7,11,12,15,17,18,10] == subjectNumber)
    sessions = 11;
end

% Area groups (the vertices of the listed areas are pooled)
groupNames = {'PMd',        'M1'};
groupAreas = {{'6a','6d'},  {'4'}};

hemispheres = {'r','l'};
nTime       = 201;

% Sessions in the output file name, e.g. '11_2_3_4_5_6'
sessionString = strjoin(arrayfun(@num2str, sessions, 'UniformOutput', false), '_');


%% =======================================================================
% LOOP OVER AREA GROUPS
%% =======================================================================

for group = 1:numel(groupNames)

    group_areas = groupAreas{group};

    %% ===================================================================
    % LOAD THE MEG DATA OF THE AREAS (both hemispheres) AND BEHAVIOR
    %
    % meg: trials x vertices x time, baseline-corrected
    %% ===================================================================

    meg          = [];
    behavior     = [];
    session_nums = [];

    for session = sessions(:)'
        for recordingNumber = Rec(:)'

            recordingMEG = [];

            for area = 1:numel(group_areas)
                for hemisphere = 1:2
                    erfFile = fullfile(erfFolder, [subject 'SES' num2str(session) '_' num2str(recordingNumber) ...
                        '_' group_areas{area} '_ERF_' hemispheres{hemisphere} '.mat']);
                    if isfile(erfFile)
                        ERF = load(erfFile);
                        if ~isempty(ERF.comb_dict{1,1}.erfdata)
                            recordingMEG = cat(2, recordingMEG, ERF.comb_dict{1,1}.erfdata);
                        end
                    end
                end
            end

            trialsFile = fullfile(trialsFolder, [subject '_' num2str(session) '_' num2str(recordingNumber) '_trials.mat']);

            if isfile(trialsFile)
                load(trialsFile);
                meg          = cat(1, meg, recordingMEG - mean(recordingMEG(:,:,1:201),3));
                behavior     = [behavior; data.trialinfo];
                session_nums = [session_nums; repmat(session, size(data.trialinfo,1), 1)];
            end

        end
    end

    % Downsample in time (factor 4)
    megDownsampled = [];
    for vertex = 1:size(meg,2)
        megDownsampled(:,vertex,:) = downsample(squeeze(meg(:,vertex,:))', 4)';
    end
    meg = megDownsampled;
    clear megDownsampled


    %% ===================================================================
    % LABELS
    %
    % trialinfo columns: 2 = choice (1 = right), 3 = confidence (1/2),
    % 4 = stimulus, 5 = RT, 16:25 = left contrasts, 26:35 = right
    %% ===================================================================

    % Confidence codes are inverted for right choices in the trial files
    rightLow  = behavior(:,2) == 1 & behavior(:,3) == 1;
    rightHigh = behavior(:,2) == 1 & behavior(:,3) == 2;
    behavior(rightLow,3)  = 2;
    behavior(rightHigh,3) = 1;

    % Choice: +1 = right, -1 = left
    Y_all = double(behavior(:,2));
    Y_all(Y_all == 0) = -1;

    evidence_strength = abs(mean(behavior(:,26:35),2) - mean(behavior(:,16:25),2));
    confidence_all    = behavior(:,3);
    rts_all           = behavior(:,5);
    stimulus_all      = behavior(:,4);


    %% ===================================================================
    % TRAINING: ONE DECODER PER TRAINING TIME POINT AND FOLD
    %
    % W(:,time,fold), B(time,fold): decoder in vertex space,
    % log odds = X * W + B
    %% ===================================================================

    nTrials = size(meg,1);
    nVert   = size(meg,2);
    folds   = cvpartition(nTrials, 'KFold', 10);
    nFolds  = folds.NumTestSets;

    % First 201 time points (-0.5 to 1.5 s from stimulus onset); z-score
    % each vertex across trials, at each time point (= normalize() of
    % decode_choice_bilateral.m)
    meg = meg(:,:,1:nTime);
    meg = (meg - mean(meg,1)) ./ std(meg,0,1);

    W = nan(nVert, nTime, nFolds);
    B = nan(nTime, nFolds);

    id_pca_all    = nan(nTime,2);
    explained_all = cell(nTime,1);

    for time = 1:nTime

        % PCA on all trials
        [coefficients, scores, ~, ~, explained, pcaMean] = pca(meg(:,:,time));

        nPCs = find(cumsum(explained) > 90, 1);
        id_pca_all(time,:)  = [nPCs size(scores,2)];
        explained_all{time} = explained;
        if isempty(nPCs), nPCs = 2; end
        nPCs = max(2, nPCs);

        pcWeights = coefficients(:,1:nPCs);     % vertices x PCs
        PCs       = scores(:,1:nPCs);           % = (X - pcaMean) * pcWeights

        for fold = 1:nFolds

            trainTrials = training(folds,fold);

            X_train = PCs(trainTrials,:);
            Y_train = Y_all(trainTrials);

            % Upsample the smaller class
            negative = find(Y_train == -1);
            positive = find(Y_train ==  1);

            if isempty(negative) || isempty(positive)
                continue
            end

            if numel(negative) > numel(positive)
                extra      = positive(randi(numel(positive), numel(negative) - numel(positive), 1));
                X_balanced = [X_train(negative,:); X_train(positive,:); X_train(extra,:)];
                Y_balanced = [-ones(numel(negative),1); ones(numel(positive),1); ones(numel(extra),1)];
            elseif numel(positive) > numel(negative)
                extra      = negative(randi(numel(negative), numel(positive) - numel(negative), 1));
                X_balanced = [X_train(negative,:); X_train(positive,:); X_train(extra,:)];
                Y_balanced = [-ones(numel(negative),1); ones(numel(positive),1); -ones(numel(extra),1)];
            else
                X_balanced = X_train;
                Y_balanced = Y_train;
            end

            % z-score with the training mean and SD, logistic regression
            [X_balanced, mu, sigma] = zscore(X_balanced);
            model = fitclinear(X_balanced, Y_balanced, 'Learner', 'logistic');

            % PCA projection + z-scoring + logistic regression in vertex space:
            % log odds = (((X - pcaMean)*pcWeights - mu) ./ sigma) * Beta + Bias
            %          = X * [pcWeights * (Beta./sigma')]
            %            + [Bias - (pcaMean*pcWeights + mu) * (Beta./sigma')]
            scaledBeta  = model.Beta ./ sigma';
            W(:,time,fold) = pcWeights * scaledBeta;
            B(time,fold)   = model.Bias - (pcaMean * pcWeights + mu) * scaledBeta;

        end

    end


    %% ===================================================================
    % TESTING: EACH DECODER AT EVERY TIME POINT OF THE HELD-OUT TRIALS
    %% ===================================================================

    meg = permute(meg, [1 3 2]);        % trials x time x vertices

    accuracy_gen_folds = nan(nFolds, nTime, nTime);   % fold x training time x testing time
    log_odds_all       = nan(nTrials, nTime);         % diagonal (training time = testing time)

    for fold = 1:nFolds

        testTrials = find(test(folds,fold));
        nTest      = numel(testTrials);

        % (trial, testing time) x vertices -> log odds of every training time
        X_test  = reshape(meg(testTrials,:,:), nTest*nTime, nVert);
        logOdds = X_test * W(:,:,fold) + B(:,fold)';
        logOdds = reshape(logOdds, nTest, nTime, nTime);      % trials x testing time x training time

        % Predicted choice (log odds = 0 counted as left, as predict())
        predictedChoice = sign(logOdds);
        predictedChoice(predictedChoice == 0) = -1;

        accuracy = squeeze(mean(predictedChoice == Y_all(testTrials), 1))';   % training time x testing time
        accuracy(isnan(B(:,fold)), :) = NaN;                                    % folds skipped in training
        accuracy_gen_folds(fold,:,:) = accuracy;

        for time = 1:nTime
            log_odds_all(testTrials,time) = logOdds(:,time,time);
        end

    end

    clear meg X_test logOdds predictedChoice

    % Mean over folds: rows = training time, columns = testing time
    accuracy_gen = squeeze(mean(accuracy_gen_folds, 1, 'omitnan'));

    % Diagonal of each fold (fold x time, as accuracy_all of decode_choice_bilateral.m)
    accuracy_all = nan(nFolds, nTime);
    for time = 1:nTime
        accuracy_all(:,time) = accuracy_gen_folds(:,time,time);
    end


    %% ===================================================================
    % SAVE
    %% ===================================================================

    if ~exist(outputFolder,'dir')
        mkdir(outputFolder);
    end

    save(fullfile(outputFolder, [subject '_' sessionString '_' groupNames{group} '_CHOICEdecoder_GENERALIZATION.mat']), ...
        'accuracy_gen','accuracy_gen_folds','accuracy_all','log_odds_all', ...
        'Y_all','stimulus_all','confidence_all','rts_all','evidence_strength', ...
        'session_nums','id_pca_all','explained_all','group_areas','nVert','-v7.3');

end

end
