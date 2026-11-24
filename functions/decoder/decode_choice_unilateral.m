function decode_choice_unilateral(subject, sessions, Rec, areas, areaLabel)
%% =======================================================================
% Single-hemisphere choice decoder (Figures 2F-I, 3, S6D, S6F, S9, S10)
%

%
% Inputs:
%   subject   - subject number (e.g. 1 for S01)
%   sessions  - session number(s)
%   Rec       - recording number(s) of the session
%   areas     - areas decoded together, e.g. {'6a'}, {'6a','6d'}, {'4'}
%   areaLabel - area name in the output file, e.g. '6a', '6A6D', '4'
%

%% =======================================================================


%% =======================================================================
% PATHS
%% =======================================================================

erfFolder    = '~/Desktop/2AFC/ERF';
trialsFolder = '~/Desktop/2AFC/trials';

functionFolder = fileparts(mfilename('fullpath'));
outputFolder   = fullfile(fileparts(fileparts(functionFolder)), 'data', 'decoder_output', 'PMd_M1_hemi');


%% =======================================================================
% SETTINGS
%% =======================================================================




subjectNumber = subject;
subject = sprintf('S%02d', subjectNumber);

if sessions == 1 && any([1,2,4,5,7,11,12,15,17,18,10] == subjectNumber)
    sessions = 11;
end

hemispheres = {'l','r'};
nTime       = 201;


%% =======================================================================
% LOOP OVER HEMISPHERES AND SESSIONS
%% =======================================================================

for hemisphere = 1:2
    for session = sessions(:)'

        %% ===============================================================
        % LOAD THE MEG DATA OF THE AREAS (this hemisphere) AND BEHAVIOR
        %
        % meg: trials x vertices x time
        %% ===============================================================

        meg      = [];
        behavior = [];

        for recordingNumber = Rec(:)'

            recordingMEG = [];

            for area = 1:numel(areas)
                erfFile = fullfile(erfFolder, [subject 'SES' num2str(session) '_' num2str(recordingNumber) ...
                    '_' areas{area} '_ERF_' hemispheres{hemisphere} '.mat']);
                if isfile(erfFile)
                    ERF = load(erfFile);
                    recordingMEG = cat(2, recordingMEG, ERF.comb_dict{1,1}.erfdata);
                end
            end

            trialsFile = fullfile(trialsFolder, [subject '_' num2str(session) '_' num2str(recordingNumber) '_trials.mat']);

            if isfile(trialsFile) && ~isempty(recordingMEG)
                load(trialsFile);
    meg      = cat(1, meg, (recordingMEG - mean(recordingMEG(:,:,1:201),3)));
                behavior = [behavior; data.trialinfo];
            end

        end

        if isempty(meg)
            continue
        end

        % Downsample in time (factor 4)
        megDownsampled = [];
        for vertex = 1:size(meg,2)
            megDownsampled(:,vertex,:) = downsample(squeeze(meg(:,vertex,:))', 4)';
        end
        meg = megDownsampled;


        %% ===============================================================
        % LABELS
        %
        % trialinfo columns: 2 = choice (1 = right), 3 = confidence (1/2),
        % 4 = stimulus, 5 = RT, 16:25 = left contrasts, 26:35 = right
        %% ===============================================================

        % Confidence codes are inverted for right choices in the trial files
        rightLow  = behavior(:,2) == 1 & behavior(:,3) == 1;
        rightHigh = behavior(:,2) == 1 & behavior(:,3) == 2;
        behavior(rightLow,3)  = 2;
        behavior(rightHigh,3) = 1;

        % Contralateral choice: +1 = right (left hemisphere) or left (right hemisphere)
        choice = double(behavior(:,2));
        choice(choice == 0) = -1;

        if hemisphere == 1
            Y_all = choice;
        else
            Y_all = -choice;
        end

        evidence_strength = abs(mean(behavior(:,26:35),2) - mean(behavior(:,16:25),2));
        confidence_all    = behavior(:,3);
        rts_all           = behavior(:,5);
        stimulus_all      = behavior(:,4);
        session_nums      = repmat(session, size(Y_all));


        %% ===============================================================
        % DECODING, TIME POINT BY TIME POINT
        %% ===============================================================

        nTrials = size(meg,1);
        folds   = cvpartition(nTrials, 'KFold', 10);

        log_odds_all  = nan(nTrials,nTime);
        accuracy_all  = nan(folds.NumTestSets,nTime);
        id_pca_all    = nan(nTime,1);
        explained_all = cell(nTime,1);

        for time = 1:nTime

            % PCA on the z-scored vertices (all trials of the session)
            [~, scores, ~, ~, explained] = pca(normalize(meg(:,:,time)));

            nPCs = find(cumsum(explained) > 90, 1);
            id_pca_all(time)    = nPCs;
            explained_all{time} = explained;
            if isempty(nPCs), nPCs = 2; end

            PCs = scores(:,1:max(2,nPCs));

            for fold = 1:folds.NumTestSets

                trainTrials = training(folds,fold);
                testTrials  = test(folds,fold);

                X_train = PCs(trainTrials,:);
                Y_train = Y_all(trainTrials);

                % Upsample the smaller class
                negative = find(Y_train == -1);
                positive = find(Y_train ==  1);

                if isempty(negative) || isempty(positive)
                    continue
                end

                if numel(negative) > numel(positive)
                    extra   = positive(randi(numel(positive), numel(negative) - numel(positive), 1));
                    Y_extra = ones(numel(extra),1);
                else
                    extra   = negative(randi(numel(negative), numel(positive) - numel(negative), 1));
                    Y_extra = -ones(numel(extra),1);
                end

                X_balanced = [X_train(negative,:); X_train(positive,:); X_train(extra,:)];
                Y_balanced = [-ones(numel(negative),1); ones(numel(positive),1); Y_extra];

                % z-score with the training mean and SD, logistic regression
                [X_balanced, mu, sigma] = zscore(X_balanced);
                model = fitclinear(X_balanced, Y_balanced, 'Learner', 'logistic');

                X_test = (PCs(testTrials,:) - mu) ./ sigma;

                log_odds_all(testTrials,time) = X_test * model.Beta + model.Bias;
                accuracy_all(fold,time)       = mean(predict(model, X_test) == Y_all(testTrials));

            end

        end


        %% ===============================================================
        % SAVE
        %% ===============================================================

        if ~exist(outputFolder,'dir')
            mkdir(outputFolder);
        end

        save(fullfile(outputFolder, [subject '_SES' num2str(session) '_' hemispheres{hemisphere} ...
            '_CHOICEdecoder_' areaLabel '_all.mat']), ...
            'log_odds_all','Y_all','stimulus_all','confidence_all','rts_all', ...
            'evidence_strength','accuracy_all','session_nums','id_pca_all','explained_all','-v7.3');

    end
end

end
