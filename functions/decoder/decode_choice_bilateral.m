function decode_choice_bilateral(subject, sessions, Rec)
%% =======================================================================
% Bilateral choice decoder of each of the 180 areas (Figures 4, 5)

% Inputs:
%   subject  - subject number (e.g. 1 for S01)
%   sessions - session number(s)
%   Rec      - recording number(s) of the session
%
% Output (one file per area):
%   data/decoder_output/ROI_bilateral/
%     <subject>_<sessions>_<area>_CHOICEdecoder_alltrials_all_b.mat
%
% Requires the source-reconstructed single-trial MEG data (ERF) and the
% trial information files (not part of the repository, see PATHS)
%% =======================================================================


%% =======================================================================
% PATHS
%% =======================================================================

erfFolder    = '~/Desktop/2AFC/ERF';
trialsFolder = '~/Desktop/2AFC/trials';

functionFolder = fileparts(mfilename('fullpath'));
outputFolder   = fullfile(fileparts(fileparts(functionFolder)), 'data', 'decoder_output', 'ROI_bilateral');


%% =======================================================================



% Subject ID; participants with an additional session have session 1
% stored as session 11
subjectNumber = subject;
subject = sprintf('S%02d', subjectNumber);

if sessions == 1 && any([1,2,4,5,7,11,12,15,17,18,10] == subjectNumber)
    sessions = 11;
end

areas = {'V1','V2','V3','V4','V6','V3A','V7','IPS1','V3B','V6A','V8','FFC','PIT','VMV1','VMV3','VMV2','VVC',...
    'MST','LO1','LO2','MT','PH','V4t','FST','V3CD','LO3','4','3b','1','2','3a','5m','5mv',...
    '23c','5L','24dd','24dv','SCEF','6ma','6mp','FEF','PEF','55b','6d','6v','6r','6a',...
    '43','OP4','OP1','OP2-3','FOP1','A1','52','RI','PFcm','PBelt','MBelt','LBelt','TA2','STGa',...
    'A5','STSda','STSdp','STSvp','A4','STSva','PoI2','FOP4','MI','Pir','AVI','AAIC','FOP3','FOP2','PoI1',...
    'Ig','FOP5','PI','EC','PreS','H','PeEc','PHA1','PHA3','TF','PHA2','TGd','TE1a','TE1p','TE2a','TE2p',...
    'PHT','TGv','TE1m','PSL','STV','TPOJ1','TPOJ2','TPOJ3','7Pm','7AL','7Am','7PL','7PC','LIPv','VIP',...
    'MIP','LIPd','AIP','PFt','PGp','IP2','IP1','IP0','PFop','PF','PFm','PGi','PGs','RSC','POS2','PCV','7m',...
    'POS1','23d','v23ab','d23ab','31pv','ProS','DVT','31pd','31a','p24pr','33pr','a24pr','p32pr','a24',...
    'd32','8BM','p32','10r','9m','10v','25','s32','pOFC','a32pr','p24','47m','10d','a10p','10pp','11l','13l',...
    'OFC','47s','p10p','44','45','47l','a47r','IFJa','IFJp','IFSp','IFSa','p47r','SFL','8Av','8Ad','8BL',...
    '9p','8C','p9-46v','46','a9-46v','9-46d','9a','i6-8','s6-8'};

hemispheres = {'r','l'};
nTime       = 201;

% Sessions in the output file name, e.g. '11_2_3_4_5_6'
sessionString = strjoin(arrayfun(@num2str, sessions, 'UniformOutput', false), '_');


%% =======================================================================
% LOOP OVER AREAS
%% =======================================================================

for area = 1:numel(areas)

    %% ===================================================================
    % LOAD THE MEG DATA OF THE AREA (both hemispheres) AND BEHAVIOR
    %
    % meg: trials x vertices x time, baseline-corrected
    %% ===================================================================

    meg          = [];
    behavior     = [];
    session_nums = [];

    for session = sessions(:)'
        for recordingNumber = Rec(:)'

            recordingMEG = [];

            for hemisphere = 1:2
                erfFile = fullfile(erfFolder, [subject 'SES' num2str(session) '_' num2str(recordingNumber) ...
                    '_' areas{area} '_ERF_' hemispheres{hemisphere} '.mat']);
                if isfile(erfFile)
                    ERF = load(erfFile);
                    recordingMEG = cat(2, recordingMEG, ERF.comb_dict{1,1}.erfdata);
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
    % DECODING, TIME POINT BY TIME POINT
    %% ===================================================================

    nTrials = size(meg,1);
    folds   = cvpartition(nTrials, 'KFold', 10);

    log_odds_all  = nan(nTrials,nTime);
    accuracy_all  = nan(folds.NumTestSets,nTime);
    id_pca_all    = nan(nTime,2);
    explained_all = cell(nTime,1);

    for time = 1:nTime

        % PCA on the z-scored vertices (all trials)
        [~, scores, ~, ~, explained] = pca(normalize(meg(:,:,time)));

        nPCs = find(cumsum(explained) > 90, 1);
        id_pca_all(time,:)  = [nPCs size(scores,2)];
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


    %% ===================================================================
    % SAVE
    %% ===================================================================

    if ~exist(outputFolder,'dir')
        mkdir(outputFolder);
    end

    save(fullfile(outputFolder, [subject '_' sessionString '_' areas{area} '_CHOICEdecoder_alltrials_all_b.mat']), ...
        'log_odds_all','Y_all','stimulus_all','confidence_all','rts_all', ...
        'evidence_strength','accuracy_all','session_nums','id_pca_all','explained_all','-v7.3');

end

end
