function Decoder_Performance = decode_choice_unilateral(subject, sessions, Rec)

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

% === Fixed aggregated areas =============================================
areas_to_decode = {'6a','6d','4'};

% Build session string for naming
sess_str = strjoin(arrayfun(@num2str, sessions, 'UniformOutput', false), '_');

% Add required paths 
addpath('~/Desktop/2AFC/Analyses/liblinear-2.45/matlab');
addpath('~/aletoso/Desktop/2AFC/trials/');

% === Hemisphere decoding ==================================================
hem_sides = {'l','r'};  % Decode left, then right separately
nH = numel(hem_sides);

for h = 1:nH
    hemi = hem_sides{h};
    fprintf('\n ====================\n Decoding hemisphere: %s\n ====================\n', hemi);

    % === Loop over sessions (SEPARATE DECODER EACH SESSION) ===============
    for si = 1:numel(sessions)
        this_sess = sessions(si);
        fprintf('\n--- SESSION %d ---\n', this_sess);

        % Reset storage for this session
        All_ERF_R = [];
        All_Behavior = [];
        session_nums = [];

        % === Load ERF + Behavior (aggregate areas) ========================
        for ri = 1:numel(Rec)
            All_ERF = [];

            for a_n = 1:numel(areas_to_decode)
                area_name = areas_to_decode{a_n};

                path = fullfile('~/Desktop/2AFC/ERF',...
                    [subject, 'SES', num2str(this_sess), '_', num2str(Rec(ri)), '_', area_name, '_ERF_', hemi, '.mat']);

                if isfile(path)
                    tmp = load(path);
                    erf = tmp.comb_dict{1,1}.erfdata;
                    if ~isempty(erf)
                        All_ERF = cat(2, All_ERF, erf);
                    end
                end
            end

            behav_path = fullfile('~/Desktop/2AFC/trials',...
                [subject, '_', num2str(this_sess), '_', num2str(Rec(ri)), '_trials.mat']);

            if isfile(behav_path)
                load(behav_path);
                All_ERF = permute(All_ERF, [3,1,2]);
                All_ERF = permute(All_ERF, [2,3,1]);
                baseline_mean = mean(All_ERF(:,:,1:201), 3);
                All_ERF_R = cat(1, All_ERF_R, All_ERF );
                All_Behavior = [All_Behavior; data.trialinfo];
                session_nums = [session_nums; repmat(this_sess, size(data.trialinfo,1), 1)];
            end
        end

        if isempty(All_ERF_R)
            fprintf('No trials in session %d\n', this_sess);
            continue
        end

        % === Downsample time ==============================================
        All_ERF_R_ds = [];
        for vi = 1:size(All_ERF_R,2)
            tmp = permute(All_ERF_R(:,vi,:), [3,1,2]);
            tmp_ds = downsample(tmp, 4);
            All_ERF_R_ds(:,:,vi) = tmp_ds;
        end
        All_ERF_R = permute(All_ERF_R_ds, [2,3,1]);
        clear All_ERF_R_ds tmp tmp_ds

        % === Behavior labels ==============================================
        % Fix confidence inversion
        id_1_1 = (All_Behavior(:,2)==1 & All_Behavior(:,3)==1);
        id_1_2 = (All_Behavior(:,2)==1 & All_Behavior(:,3)==2);
        All_Behavior(id_1_1,3)=2;
        All_Behavior(id_1_2,3)=1;

        % Extract choice
        choice_raw = double(All_Behavior(:,2));  % Right=1, Left=0
        choice_raw(choice_raw == 0) = -1;        % convert Left -> -1

        % === CONTRALATERAL re-labeling ====================================
        if strcmp(hemi,'l')
            % Left hemisphere decodes RIGHT choice (contralateral)
            % So RIGHT trials = 1, LEFT = -1  ==> unchanged
            Y_all = choice_raw;
        else
            % Right hemisphere decodes LEFT choice (contralateral)
            % So LEFT must be labeled as +1, RIGHT as -1
            Y_all = -choice_raw;
        end

        % Evidence strength
        left_samples  = All_Behavior(:,16:25);
        right_samples = All_Behavior(:,26:35);
        evidence_strength = abs(mean(right_samples,2) - mean(left_samples,2));

        % Save behavioral vars
        confidence_all = All_Behavior(:,3);
        rts_all        = All_Behavior(:,5);
        stimulus_all   = All_Behavior(:,4);

        % === Cross-validation ===============
        nTrials = size(All_ERF_R,1);
        Time   = 1:201;
        nTime  = numel(Time);

        cv = cvpartition(nTrials, 'KFold', 10);
        log_odds_all = nan(nTrials, nTime);
        accuracy_all = nan(cv.NumTestSets, nTime);
id_pca_all=nan(nTime,1);
        for t = 1:nTime
            X = squeeze(All_ERF_R(:,:,Time(t)));
            X = normalize(X);
% PCA, get N PCA that explains 90% of the variance
            [~, score, ~, ~, explained] = pca(X);
            cumexp = cumsum(explained);
                    explained_all{t,:}=explained;

            id_pca = find(cumexp > 90, 1);
            id_pca_all(t,:)=id_pca;
            if isempty(id_pca), id_pca = 2; end
            X_all = score(:,1:max(2,id_pca));

            for fold = 1:cv.NumTestSets
                train_idx = training(cv, fold);
                test_idx  = test(cv, fold);

                Y_train = Y_all(train_idx);
                X_train = X_all(train_idx,:);

              % Identify classes
idL = find(Y_train == -1);
idR = find(Y_train == 1);

nL = numel(idL);
nR = numel(idR);

% If one class is empty, skip fold
if nL < 1 || nR < 1
    continue;
end


% % --- UPSAMPLING IMPLEMENTATION ---
if nL > nR
    % Right is minority . upsample Right
    n_extra = nL - nR;
    id_extra = idR(randi(nR, n_extra, 1));

    X_train_bal = [X_train(idL,:); X_train(idR,:); X_train(id_extra,:)];
    Y_train_bal = [-1*ones(nL,1);  ones(nR,1);    ones(n_extra,1)];

else
    % Left is minority - upsample Left
    n_extra = nR - nL;
    id_extra = idL(randi(nL, n_extra, 1));

    X_train_bal = [X_train(idL,:); X_train(idR,:); X_train(id_extra,:)];
    Y_train_bal = [-1*ones(nL,1);  ones(nR,1);   -1*ones(n_extra,1)];
end

                [X_train_z, mu, sigma] = zscore(X_train_bal);
                Mdl = fitclinear(X_train_z, Y_train_bal, 'Learner','logistic');

                X_test_z = (X_all(test_idx,:)-mu)./sigma;
                Y_pred   = predict(Mdl, X_test_z);
                logodds  = X_test_z*Mdl.Beta + Mdl.Bias;

                log_odds_all(test_idx, t) = logodds;
                accuracy_all(fold, t)    = mean(Y_pred == Y_all(test_idx));
            end
        end

        % === SAVE OUTPUT SEPARATELY FOR EACH HEMISPHERE & SESSION ========
        savepath = '/home/aletoso/Desktop/2AFC/Figures_confidence/wnag/NO_B/';
        if ~exist(savepath,'dir'), mkdir(savepath); end

        filename = fullfile(savepath,...
            [subject, '_SES', num2str(this_sess), '_', hemi, '_CHOICEdecoder_6A6D4_all_nob_80.mat']);

        save(filename, 'log_odds_all','Y_all','stimulus_all','confidence_all',...
            'rts_all','evidence_strength','accuracy_all','session_nums','id_pca_all','explained_all','-v7.3');

    end
end

Decoder_Performance = [];
end
