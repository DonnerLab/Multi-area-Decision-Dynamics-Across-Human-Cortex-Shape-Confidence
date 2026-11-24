function decode_choice_all_betas(subject, sessions, Rec)
%% =======================================================================
% Choice decoder on all vertices of PMd/M1 (areas 4, 6a, 6d, both
% hemispheres): decoder weights of each vertex (Figure 2E)
%
% Inputs:
%   subject  - subject number (e.g. 1 for S01)
%   sessions - session number
%   Rec      - recording number(s) of the session
%
% Output:
%   data/decoder_output/all_VX/
%     <subject>_train_CHOICE_decoder_session<n>_PCA_PMd_M1_SINGLETRIAL.mat
%     (Prob_Performance: 1 x time x vertices, decoder weights)
%
%% =======================================================================


%% =======================================================================
% PATHS
%% =======================================================================

erfFolder    = '~/Desktop/2AFC/ERF';
trialsFolder = '~/Desktop/2AFC/trials';

functionFolder = fileparts(mfilename('fullpath'));
outputFolder   = fullfile(fileparts(fileparts(functionFolder)), 'data', 'decoder_output', 'all_VX');


%% =======================================================================
% SETTINGS
%% =======================================================================

% Subject ID; participants with an additional session have session 1
% stored as session 11
subjectNumber = subject;
subject = sprintf('S%02d', subjectNumber);

if sessions == 1 && any([1,2,4,5,7,11,12,15,17,18,10] == subjectNumber)
    sessions = 11;
end

areas       = {'4','6a','6d'};
hemispheres = {'l','r'};


%% =======================================================================
% LOAD THE MEG DATA AND THE BEHAVIOR
%
% meg: trials x vertices x time
%% =======================================================================

meg      = [];
behavior = [];

for recordingNumber = Rec(:)'

    if ~isfile(fullfile(erfFolder, [subject 'SES' num2str(sessions) '_' num2str(recordingNumber) '_' areas{1} '_ERF_r.mat']))
        continue
    end

    recordingMEG = [];

    for hemisphere = 1:2

        for area = 1:numel(areas)
            erfFile = fullfile(erfFolder, [subject 'SES' num2str(sessions) '_' num2str(recordingNumber) ...
                '_' areas{area} '_ERF_' hemispheres{hemisphere} '.mat']);
            if isfile(erfFile)
                ERF = load(erfFile);
                recordingMEG = cat(2, ERF.comb_dict{1,1}.erfdata, recordingMEG);
            end
        end

        % Number of vertices after this hemisphere
        n_vx_side(hemisphere,:) = size(recordingMEG,2);

    end

    load(fullfile(trialsFolder, [subject '_' num2str(sessions) '_' num2str(recordingNumber) '_trials.mat']));

    %  baseline-corrected activity
    meg      = cat(1, meg, (recordingMEG - mean(recordingMEG(:,:,1:201),3)));
    behavior = [behavior; data.trialinfo];

end


%% =======================================================================
% DECODING, TIME POINT BY TIME POINT
%% =======================================================================

% Choice: +1 = right, -1 = left
Y_all = double(behavior(:,2));
Y_all(Y_all == 0) = -1;

folds = cvpartition(size(behavior,1), 'Kfold', 10);

nTime = size(meg,3);

Prob_Performance    = nan(1,nTime,size(meg,2));   % decoder weights
Decoder_Performance = nan(1,nTime);               % decoding accuracy

for time = 1:nTime

    X_all = meg(:,:,time);

    predictedChoice = [];
    testChoice      = [];
    foldWeights     = [];

    for fold = 1:10

        trainTrials = training(folds,fold);
        testTrials  = test(folds,fold);

        X = X_all(trainTrials,:);
        Y = Y_all(trainTrials);

        % Upsample the smaller class
        negative = find(Y == -1);
        positive = find(Y ==  1);

        if numel(negative) > numel(positive)
            extra = positive(randi(numel(positive), numel(negative) - numel(positive), 1));
        else
            extra = negative(randi(numel(negative), numel(positive) - numel(negative), 1));
        end

        X = [X; X(extra,:)];
        Y = [Y; Y(extra)];

        % z-score with the training mean and SD, logistic regression
        [X, mu, sigma] = zscore(X);
        model = fitclinear(X, Y, 'Learner', 'logistic');

        predictedChoice = [predictedChoice; predict(model, normalize(X_all(testTrials,:), 'center', mu, 'scale', sigma))];
        testChoice      = [testChoice; Y_all(testTrials)];
        foldWeights     = [foldWeights; model.Beta'];

    end

    Decoder_Performance(time)  = mean(predictedChoice == testChoice);
    Prob_Performance(1,time,:) = mean(foldWeights,1);

end


%% =======================================================================
% SAVE
%% =======================================================================

if ~exist(outputFolder,'dir')
    mkdir(outputFolder);
end

All_samples  = behavior(:,16:35);
All_Choices  = behavior(:,2);
block_type   = behavior(:,13);
All_Behavior = behavior;

save(fullfile(outputFolder, [subject '_train_CHOICE_decoder_session' num2str(sessions) '_PCA_PMd_M1_SINGLETRIAL.mat']), ...
    'Prob_Performance','Decoder_Performance','All_samples','All_Choices','block_type','All_Behavior','n_vx_side')

end
