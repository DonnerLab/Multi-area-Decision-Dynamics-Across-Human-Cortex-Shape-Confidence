function Decoder_Performance = decode_choice_bilateral(subject, sessions, Rec)

% === Setup ===============================================================
rng(42);

Extra_sessions = [1,2,4,5,7,0,11,12,15,17,18,10,0,0,0,0,0,0,0,0,0];
s_n = subject;

if s_n < 10
    subject = ['S0', num2str(s_n)];
else
    subject = ['S', num2str(s_n)];
end

if sessions == 1 && any(Extra_sessions == s_n)
    sessions = 11;
end

% Areas to decode
areas_to_decode = {'V1','V2','V3','V4','V6','V3A','V7','IPS1','V3B','V6A','V8','FFC','PIT','VMV1','VMV3','VMV2','VVC',...
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

% Build session string for saving
sess_str = strjoin(arrayfun(@num2str, sessions, 'UniformOutput', false), '_');

% Add required paths
addpath('~/Desktop/2AFC/Analyses/liblinear-2.45/matlab');
addpath('~/aletoso/Desktop/2AFC/trials/');

% === Loop over areas =====================================================
for a_n = 1:numel(areas_to_decode)
    area_name = areas_to_decode{a_n};
    fprintf('\n>>> Decoding CHOICE for %s - %s <<<\n', subject, area_name);

    All_ERF_R = [];
    All_Behavior = [];
    session_nums = [];

    sides = {'r','l'};
   

    % === Load MEG data + Behavior =============================================
    for si = 1:numel(sessions)
        for ri = 1:numel(Rec)
            All_ERF = [];
            for side_i = 1:2
                path = fullfile('~/Desktop/2AFC/ERF', ...
                    [subject, 'SES', num2str(sessions(si)), '_', num2str(Rec(ri)), '_', area_name, '_ERF_', sides{side_i}, '.mat']);
                if isfile(path)
                    tmp = load(path);
                    erf = tmp.comb_dict{1,1}.erfdata;
                    if ~isempty(erf)
                        All_ERF = cat(2, All_ERF, erf);
                    end
                end
            end

            behav_path = fullfile('~/Desktop/2AFC/trials', ...
                [subject, '_', num2str(sessions(si)), '_', num2str(Rec(ri)), '_trials.mat']);
            if isfile(behav_path)
                load(behav_path);
                All_ERF = permute(All_ERF, [3,1,2]);
                All_ERF = permute(All_ERF, [2,3,1]);
                % baseine subtracting pre-stimulus 500 ms per each trial
                baseline_mean = mean(All_ERF(:,:,1:201), 3);
                All_ERF_R = cat(1, All_ERF_R, All_ERF - baseline_mean );
                All_Behavior = [All_Behavior; data.trialinfo];
                session_nums = [session_nums; repmat(sessions(si), size(data.trialinfo,1), 1)];
            end
        end
    end

    % === Downsample time dimension =======================================
    All_ERF_R_ds = [];
    for vi = 1:size(All_ERF_R,2)
        tmp = permute(All_ERF_R(:,vi,:), [3,1,2]);
        tmp_ds = downsample(tmp, 4);
        All_ERF_R_ds(:,:,vi) = tmp_ds;
    end
    All_ERF_R = permute(All_ERF_R_ds, [2,3,1]);
    clear All_ERF_R_ds tmp tmp_ds

    % === Prepare behavioral labels =======================================
    % Fix confidence inversion
    id_1_1 = (All_Behavior(:,2)==1 & All_Behavior(:,3)==1);
    id_1_2 = (All_Behavior(:,2)==1 & All_Behavior(:,3)==2);
    All_Behavior(id_1_1,3)=2;
    All_Behavior(id_1_2,3)=1;

    % Choice labels (to decode)
    Y_all = double(All_Behavior(:,2)); % choice column
    Y_all(Y_all == 0) = -1;

    % Evidence strength
    left_samples  = All_Behavior(:,16:25);
    right_samples = All_Behavior(:,26:35);
    evidence_strength = abs(mean(right_samples,2) - mean(left_samples,2));

    % Store other variables for saving
    confidence_all = All_Behavior(:,3);
    rts_all        = All_Behavior(:,5);
    stimulus_all   = All_Behavior(:,4);

    % === Cross-validation setup ==========================================
    nTrials = size(All_ERF_R,1);
    Time   = 1:201;
    nTime   = size(1:201,2);
    
    cv = cvpartition(nTrials, 'KFold', 10);

    log_odds_all = nan(nTrials, nTime);
    accuracy_all = nan(cv.NumTestSets, nTime);

    % === Decode over time ================================================
    for t = 1:nTime
        X = squeeze(All_ERF_R(:,:,Time(t)));
        X = normalize(X);

        % PCA and take N PCA that eaxplain 90% of  the variance
        [~, score, ~, ~, explained] = pca(X);
        cumexp = cumsum(explained);
       
        id_pca = find(cumexp > 90, 1);
        id_pca_all(t,:)=[id_pca size(score,2)];
        explained_all{t,:}=explained;
        if isempty(id_pca), id_pca = 2; end
        X_all = score(:,1:max(2,id_pca));

        for fold = 1:cv.NumTestSets
            train_idx = training(cv, fold);
            test_idx  = test(cv, fold);

            Y_train = Y_all(train_idx);
            X_train = X_all(train_idx,:);

            % Balance training set (upsample)
           idL = find(Y_train == -1);
idR = find(Y_train == 1);

nL = numel(idL);
nR = numel(idR);

if nL < 1 || nR < 1
    continue;
end

if nL > nR
    % Right class is minority
    n_extra = nL - nR;
    id_extra = idR(randi(nR,n_extra,1));

    X_train_bal = [X_train(idL,:);
                   X_train(idR,:);
                   X_train(id_extra,:)];

    Y_train_bal = [-ones(nL,1);
                    ones(nR,1);
                    ones(n_extra,1)];

elseif nR > nL
    % Left class is minority
    n_extra = nR - nL;
    id_extra = idL(randi(nL,n_extra,1));

    X_train_bal = [X_train(idL,:);
                   X_train(idR,:);
                   X_train(id_extra,:)];

    Y_train_bal = [-ones(nL,1);
                    ones(nR,1);
                   -ones(n_extra,1)];

else
    X_train_bal = X_train;
    Y_train_bal = Y_train;
end




            [X_train_z, mu, sigma] = zscore(X_train_bal);
            Mdl = fitclinear(X_train_z, Y_train_bal, 'Learner', 'logistic');

            % Test
            X_test = X_all(test_idx,:);
            X_test_z = (X_test - mu) ./ sigma;
            Y_test = Y_all(test_idx);

            Y_pred = predict(Mdl, X_test_z);
            logodds = X_test_z * Mdl.Beta + Mdl.Bias;

            log_odds_all(test_idx, t) = logodds;
            accuracy_all(fold, t) = mean(Y_pred == Y_test);
        end
    end

    % === Save outputs ====================================================
    savepath = '/home/aletoso/Desktop/2AFC/Figures_confidence/wnag/NO_B/all_areas/';
    if ~exist(savepath, 'dir')
        mkdir(savepath);
    end

    filename = fullfile(savepath, [subject, '_', sess_str, '_', area_name, '_CHOICEdecoder_alltrials_all_b.mat']);
    save(filename, 'log_odds_all', 'Y_all', 'stimulus_all', 'confidence_all', ...
        'rts_all', 'evidence_strength', 'accuracy_all', 'session_nums','id_pca_all','explained_all','-v7.3');

end

Decoder_Performance = [];
end
