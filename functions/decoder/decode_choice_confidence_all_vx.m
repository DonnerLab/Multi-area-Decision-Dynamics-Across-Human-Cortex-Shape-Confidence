function decode_choice_confidence_all_vx(subject, sessions, Rec)
%% =======================================================================
% Choice and confidence decoders trained on all cortical areas of both
% hemispheres (Figure S4E-G)
%
%
% Labels: choice +1 = right, -1 = left; confidence +1 = high, -1 = low
%
% Inputs:
%   subject  - subject number (e.g. 1 for S01)
%   sessions - session number(s)
%   Rec      - recording number(s) of the session
%

%% =======================================================================


%% =======================================================================
% PATHS
%% =======================================================================

% Source-reconstructed single-trial MEG data and trial information
erfFolder    = '~/Desktop/2AFC/ERF';
trialsFolder = '~/Desktop/2AFC/trials';

% Output folder: <repository>/data/decoder_output/all_areas
functionFolder = fileparts(mfilename('fullpath'));
outputFolder   = fullfile(fileparts(fileparts(functionFolder)), 'data', 'decoder_output', 'all_areas');


%% =======================================================================
% SETTINGS
%% =======================================================================


% Subject ID; participants with an additional session have session 1
% stored as session 11
subjectNumber = subject;
subject = sprintf('S%02d', subjectNumber);

extraSessions = [1,2,4,5,7,11,12,15,17,18,10];

if sessions == 1 && any(extraSessions == subjectNumber)
    sessions = 11;
end

areas = { ...
    'V1','V2','V3','V4','V6','V3A','V7','IPS1','V3B','V6A','V8','FFC',...
    'PIT','VMV1','VMV3','VMV2','VVC',...
    'MST','LO1','LO2','MT','PH','V4t','FST','V3CD','LO3',...
    '4','3b','1','2','3a','5m','5mv',...
    '23c','5L','24dd','24dv','SCEF','6ma','6mp','FEF','PEF','55b',...
    '6d','6v','6r','6a',...
    '43','OP4','OP1','OP2-3','FOP1','A1','52','RI','PFcm','PBelt',...
    'MBelt','LBelt','TA2','STGa',...
    'A5','STSda','STSdp','STSvp','A4','STSva','PoI2','FOP4','MI',...
    'Pir','AVI','AAIC','FOP3','FOP2','PoI1',...
    'Ig','FOP5','PI','EC','PreS','H','PeEc','PHA1','PHA3','TF','PHA2',...
    'TGd','TE1a','TE1p','TE2a','TE2p',...
    'PHT','TGv','TE1m','PSL','STV','TPOJ1','TPOJ2','TPOJ3',...
    '7Pm','7AL','7Am','7PL','7PC','LIPv','VIP',...
    'MIP','LIPd','AIP','PFt','PGp','IP2','IP1','IP0','PFop','PF',...
    'PFm','PGi','PGs','RSC','POS2','PCV','7m',...
    'POS1','23d','v23ab','d23ab','31pv','ProS','DVT','31pd','31a',...
    'p24pr','33pr','a24pr','p32pr','a24',...
    'd32','8BM','p32','10r','9m','10v','25','s32','pOFC','a32pr',...
    'p24','47m','10d','a10p','10pp','11l','13l',...
    'OFC','47s','p10p','44','45','47l','a47r','IFJa','IFJp','IFSp',...
    'IFSa','p47r','SFL','8Av','8Ad','8BL',...
    '9p','8C','p9-46v','46','a9-46v','9-46d','9a','i6-8','s6-8'};

hemispheres = {'l','r'};

downsamplingFactor = 8;
explainedVariance  = 0.90;   % PCs explaining 90% of the training variance
nFolds             = 10;


%% =======================================================================
% LOOP OVER SESSIONS
%% =======================================================================

for session = sessions(:)'

    %% ===================================================================
    % LOAD THE MEG DATA OF ALL AREAS AND THE BEHAVIOR (all recordings)
    %
    % meg: trials x time x vertices of all areas and both hemispheres
    %% ===================================================================

    meg       = [];
    behavior  = [];
    recording = [];

    for recordingNumber = Rec(:)'

        recordingMEG = [];

        for area = 1:numel(areas)
            for hemisphere = 1:2

                erfFile = fullfile(erfFolder, [subject 'SES' num2str(session) '_' num2str(recordingNumber) ...
                    '_' areas{area} '_ERF_' hemispheres{hemisphere} '.mat']);

                if ~isfile(erfFile)
                    continue
                end

                ERF = load(erfFile);
                erf = ERF.comb_dict{1,1}.erfdata;    % trials x vertices x time

                if isempty(erf)
                    continue
                end

                % trials x time x vertices, baseline correction (first 201
                % samples = pre-stimulus 500 ms)
                erf = permute(erf, [1 3 2]);
                erf = erf - mean(erf(:,1:201,:),2);

                % Downsample in time
                erfDownsampled = zeros(size(erf,1), floor((size(erf,2)-1)/downsamplingFactor) + 1, size(erf,3));
                for vertex = 1:size(erf,3)
                    erfDownsampled(:,:,vertex) = downsample(erf(:,:,vertex)', downsamplingFactor)';
                end

                recordingMEG = cat(3, recordingMEG, erfDownsampled);

            end
        end

        if isempty(recordingMEG)
            continue
        end

        % Trial information of this recording
        load(fullfile(trialsFolder, [subject '_' num2str(session) '_' num2str(recordingNumber) '_trials.mat']));

        if size(recordingMEG,1) ~= size(data.trialinfo,1)
            error('%s session %d recording %d: MEG and behavior have different numbers of trials', ...
                subject, session, recordingNumber);
        end

        meg       = cat(1, meg, recordingMEG);
        behavior  = [behavior; data.trialinfo];
        recording = [recording; repmat(recordingNumber, size(data.trialinfo,1), 1)];

    end

    if isempty(meg)
        continue
    end


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

    % Confidence: +1 = high (2), -1 = low (1); other codes excluded
    confidence_all = behavior(:,3);
    conf_valid     = confidence_all == 1 | confidence_all == 2;
    Conf_Y         = nan(size(confidence_all));
    Conf_Y(confidence_all == 1) = -1;
    Conf_Y(confidence_all == 2) =  1;

    evidence_strength = abs(mean(behavior(:,26:35),2) - mean(behavior(:,16:25),2));
    rts_all           = behavior(:,5);
    stimulus_all      = behavior(:,4);
    session_nums      = repmat(session, size(Y_all));
    rec_nums          = recording;

    % trials x vertices x time
    meg = permute(meg, [1 3 2]);

    [nTrials, nVertices, nTime] = size(meg);

    % Largest possible number of PCs (no upper limit is imposed)
    maxPCs = min(nVertices, nTrials);

    % Same cross-validation folds at all time points
    folds = cvpartition(nTrials, 'KFold', nFolds);


    %% ===================================================================
    % OUTPUTS
    %% ===================================================================

    % Cross-validated log odds of each trial
    log_odds_choice = nan(nTrials,nTime);
    log_odds_conf   = nan(nTrials,nTime);

    % Accuracy of each fold
    accuracy_choice     = nan(nFolds,nTime);
    bal_accuracy_choice = nan(nFolds,nTime);
    accuracy_conf       = nan(nFolds,nTime);
    bal_accuracy_conf   = nan(nFolds,nTime);

    % Number of PCs and their explained variance
    n_pcs_all     = nan(nFolds,nTime);
    explained_all = cell(nTime,nFolds);

    % Weights in PC space (NaN after the last PC used) and in vertex space
    beta_choice_pc   = nan(maxPCs,nTime,nFolds);
    beta_conf_pc     = nan(maxPCs,nTime,nFolds);
    bias_choice      = nan(nTime,nFolds);
    bias_conf        = nan(nTime,nFolds);
    beta_choice_feat = nan(nVertices,nTime,nFolds);
    beta_conf_feat   = nan(nVertices,nTime,nFolds);

    % Similarity of the two weight vectors
    cos_sim   = nan(nFolds,nTime);
    angle_deg = nan(nFolds,nTime);

    fold_recordings = (1:nFolds)';


    %% ===================================================================
    % DECODING, TIME POINT BY TIME POINT
    %% ===================================================================

    for time = 1:nTime

        % z-score each vertex within each recording
        X = meg(:,:,time);

        for recordingNumber = unique(recording)'
            trials = recording == recordingNumber;
            sigma  = std(X(trials,:),0,1);
            sigma(sigma == 0) = 1;
            X(trials,:) = (X(trials,:) - mean(X(trials,:),1)) ./ sigma;
        end

        for fold = 1:nFolds

            trainTrials = training(folds,fold);
            testTrials  = test(folds,fold);

            % z-score with the training mean and SD
            mu    = mean(X(trainTrials,:),1);
            sigma = std(X(trainTrials,:),0,1);
            sigma(sigma == 0) = 1;

            X_train = (X(trainTrials,:) - mu) ./ sigma;
            X_test  = (X(testTrials,:)  - mu) ./ sigma;

            % PCA on the training trials: PCs explaining 90% of the variance
            [coefficients,~,~,~,explained] = pca(X_train, 'Centered', false);

            nPCs = find(cumsum(explained) >= 100*explainedVariance, 1);
            if isempty(nPCs), nPCs = numel(explained); end
            nPCs = min(nPCs, size(coefficients,2));

            PC_train = X_train * coefficients(:,1:nPCs);
            PC_test  = X_test  * coefficients(:,1:nPCs);

            n_pcs_all(fold,time)     = nPCs;
            explained_all{time,fold} = explained(1:nPCs);

            % ---------------------------------------------------------
            % Choice decoder
            % ---------------------------------------------------------

            [beta_c, bias_c] = balanced_logistic_regression(PC_train, Y_all(trainTrials));

            if ~isempty(beta_c)

                logOdds    = PC_test * beta_c + bias_c;
                prediction = sign(logOdds);
                prediction(prediction == 0) = -1;

                log_odds_choice(testTrials,time) = logOdds;
                accuracy_choice(fold,time)       = mean(prediction == Y_all(testTrials));
                bal_accuracy_choice(fold,time)   = balanced_accuracy(prediction, Y_all(testTrials));

                beta_choice_pc(1:nPCs,time,fold) = beta_c;
                bias_choice(time,fold)           = bias_c;
                beta_choice_feat(:,time,fold)    = coefficients(:,1:nPCs) * beta_c;

            end

            % ---------------------------------------------------------
            % Confidence decoder (same PCs; trials with a confidence
            % report)
            % ---------------------------------------------------------

            trainConfidence = conf_valid(trainTrials);
            testConfidence  = conf_valid(testTrials);
            Y_train = Conf_Y(trainTrials);
            Y_test  = Conf_Y(testTrials);

            [beta_f, bias_f] = balanced_logistic_regression(PC_train(trainConfidence,:), Y_train(trainConfidence));

            if ~isempty(beta_f)

                logOdds    = PC_test * beta_f + bias_f;
                prediction = sign(logOdds(testConfidence));
                prediction(prediction == 0) = -1;

                log_odds_conf(testTrials,time) = logOdds;
                accuracy_conf(fold,time)       = mean(prediction == Y_test(testConfidence));
                bal_accuracy_conf(fold,time)   = balanced_accuracy(prediction, Y_test(testConfidence));

                beta_conf_pc(1:nPCs,time,fold) = beta_f;
                bias_conf(time,fold)           = bias_f;
                beta_conf_feat(:,time,fold)    = coefficients(:,1:nPCs) * beta_f;

            end

            % ---------------------------------------------------------
            % Cosine similarity of the two weight vectors (identical in PC
            % space and in vertex space, since the PCs are orthonormal)
            % ---------------------------------------------------------

            if ~isempty(beta_c) && ~isempty(beta_f)
                cosine = (beta_c' * beta_f) / (norm(beta_c) * norm(beta_f));
                cos_sim(fold,time)   = max(-1, min(1, cosine));
                angle_deg(fold,time) = acosd(cos_sim(fold,time));
            end

        end

    end


    %% ===================================================================
    % SAVE
    %% ===================================================================

    if ~exist(outputFolder,'dir')
        mkdir(outputFolder);
    end

    save(fullfile(outputFolder, [subject '_SES' num2str(session) ...
        '_decoder_all_vx.mat']), ...
        'log_odds_choice','log_odds_conf', ...
        'accuracy_choice','bal_accuracy_choice', ...
        'accuracy_conf','bal_accuracy_conf', ...
        'beta_choice_pc','beta_conf_pc', ...
        'beta_choice_feat','beta_conf_feat', ...
        'bias_choice','bias_conf', ...
        'cos_sim','angle_deg', ...
        'n_pcs_all','explained_all', ...
        'fold_recordings','rec_nums', ...
        'Y_all','Conf_Y','conf_valid', ...
        'stimulus_all','confidence_all','rts_all', ...
        'evidence_strength','session_nums', ...
        '-v7.3');

end

end


function [beta, bias] = balanced_logistic_regression(X, y)
% Logistic regression (fitclinear) on a class-balanced random subsample:
% the larger class is undersampled to the size of the smaller class.
% Empty outputs if a class has fewer than 2 trials.

classNegative = find(y == -1);
classPositive = find(y ==  1);
nPerClass     = min(numel(classNegative), numel(classPositive));

if nPerClass < 2
    beta = [];
    bias = [];
    return
end

balanced = [classNegative(randperm(numel(classNegative), nPerClass)); ...
            classPositive(randperm(numel(classPositive), nPerClass))];

model = fitclinear(X(balanced,:), y(balanced), 'Learner', 'logistic');

beta = model.Beta;
bias = model.Bias;

end


function accuracy = balanced_accuracy(prediction, y)
% Mean of the accuracies of the two classes (classes present in y)

classAccuracy = [];

if any(y == -1), classAccuracy(end+1) = mean(prediction(y == -1) == -1); end
if any(y ==  1), classAccuracy(end+1) = mean(prediction(y ==  1) ==  1); end

accuracy = mean(classAccuracy);

end
