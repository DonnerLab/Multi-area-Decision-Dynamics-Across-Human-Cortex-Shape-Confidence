%% ============================================================
% PMd-projected AvgWinDV + AvgLoseDV
% Sliding-window regression
% Confidence + Reaction Time
%

%% ============================================================

clear;
close all;
clc;
rng(42);

%% ============================================================
% PATHS
%% ============================================================

scriptFolder  = fileparts(mfilename('fullpath'));
projectFolder = fileparts(fileparts(scriptFolder));

data_path = fullfile(projectFolder, ...
    'data', ...
    'decoder_output', ...
    'all_areas');

nPerm = 10000;

%% ============================================================
% SUBJECTS
%% ============================================================

subjects = { ...
    'S01','S04','S06','S07','S12','S15','S17','S18', ...
    'S11','S05','S19','S20','S22','S10','S02','S08', ...
    'S24','S21','S23','S25'};

nSubjects = numel(subjects);

%% ============================================================
% TIME SETTINGS
%% ============================================================

win_length = 35;
step_size  = 10;

win_starts = 1:step_size:(200 - win_length + 1);
nTime      = numel(win_starts);

time_sec = ((win_starts + win_length/2 - 1) * 10 - 500) / 1000;

%% ============================================================
% AREA SETS
% Only BEST 25 ROIs
%% ============================================================

areas = { ...
    'AIP','FEF','5mv','23c','i6-8','24dd','6mp','PFt','IP2', ...
    '8Av','p24pr','7PC','23d','PF','6ma','55b','7AL','24dv', ...
    'DVT','PFm','31a','RSC','5L','V6','8C'};

refAreas = {'6a','6d','4'};

set_name = 'Best 25 ROIs';

%% ============================================================
% STORAGE
%% ============================================================

betas_conf_time = cell(1, nTime);
betas_rt_time   = cell(1, nTime);

%% ============================================================
% MAIN LOOP
%% ============================================================

for si = 1:nSubjects

    subj = subjects{si};
    fprintf('Processing %s\n', subj);

    subj_data   = struct();
    valid_areas = {};
    valid_ref   = {};

    %% --------------------------------------------------------
    % LOAD ALL AREA FILES
    %% --------------------------------------------------------

    allAreas = unique([areas refAreas]);

    for ai = 1:numel(allAreas)

        area = allAreas{ai};

        files = dir(fullfile( ...
            data_path, ...
            sprintf('%s_*_%s_CHOICEdecoder_alltrials_all.mat', ...
            subj, area)));

        if isempty(files)
            continue
        end

        L  = [];
        Y  = [];
        CF = [];
        RT = [];
        EV = [];
        SS = [];

        for f = 1:numel(files)

            D = load(fullfile(files(f).folder, files(f).name));

            X = D.log_odds_all(:,1:200);

            L  = [L; X];
            Y  = [Y; D.Y_all(:)];
            CF = [CF; D.confidence_all(:) - 1];
            RT = [RT; D.rts_all(:)];
            EV = [EV; D.evidence_strength(:)];
            SS = [SS; f * ones(size(X,1),1)];

        end

        A = matlab.lang.makeValidName(area);

        subj_data.(A).logodds = L;
        subj_data.(A).choice  = Y;
        subj_data.(A).conf    = CF;
        subj_data.(A).rt      = RT;
        subj_data.(A).evid    = EV;
        subj_data.(A).sess    = SS;

        if ismember(area, areas)
            valid_areas{end+1} = A;
        end

        if ismember(area, refAreas)
            valid_ref{end+1} = A;
        end
    end

    if isempty(valid_areas) || isempty(valid_ref)
        continue
    end

    %% --------------------------------------------------------
    % COMMON SESSIONS
    %% --------------------------------------------------------

    all_names = [valid_areas valid_ref];

    sess_all = cellfun( ...
        @(a) unique(subj_data.(a).sess), ...
        all_names, ...
        'UniformOutput', false);

    common_sess = sess_all{1};

    for k = 2:numel(sess_all)
        common_sess = intersect(common_sess, sess_all{k});
    end

    if isempty(common_sess)
        continue
    end

    %% --------------------------------------------------------
    % RESTRICT TO COMMON SESSIONS
    %% --------------------------------------------------------

    for k = 1:numel(all_names)

        A = all_names{k};

        keep = ismember( ...
            subj_data.(A).sess, ...
            common_sess);

        fn = fieldnames(subj_data.(A));

        for ff = 1:numel(fn)
            subj_data.(A).(fn{ff}) = ...
                subj_data.(A).(fn{ff})(keep,:);
        end
    end

    %% --------------------------------------------------------
    % SHARED VARIABLES
    %% --------------------------------------------------------

    ref = valid_areas{1};

    choice = subj_data.(ref).choice;
    conf   = subj_data.(ref).conf;
    rt     = subj_data.(ref).rt;
    evid   = subj_data.(ref).evid;
    sess   = subj_data.(ref).sess;

    choiceSign = sign(choice);

    %% ========================================================
    % SLIDING WINDOW REGRESSION
    %% ========================================================

    for tt = 1:nTime

        idx_t = win_starts(tt):(win_starts(tt) + win_length - 1);

        %% ----------------------------------------------------
        % Compute AvgWinDV + AvgLoseDV
        %% ----------------------------------------------------

        [AvgWinDV, AvgLoseDV] = ...
            compute_AvgWin_AvgLose_timepointwise( ...
            subj_data, ...
            valid_areas, ...
            valid_ref, ...
            idx_t, ...
            choiceSign, ...
            sess);

        %% ----------------------------------------------------
        % Session-wise z-scoring
        %% ----------------------------------------------------

        AvgWin_z  = nan(size(AvgWinDV));
        AvgLose_z = nan(size(AvgLoseDV));
        evid_z    = nan(size(evid));
        rt_z      = nan(size(rt));

        for ss = unique(sess)'

            idx = sess == ss;

            valid = idx & ...
                ~isnan(AvgWinDV)  & ...
                ~isnan(AvgLoseDV) & ...
                ~isnan(evid)      & ...
                ~isnan(rt)        & ...
                ~isnan(conf);

            if sum(valid) > 2

                AvgWin_z(valid)  = zscore(AvgWinDV(valid));
                AvgLose_z(valid) = zscore(AvgLoseDV(valid));
                evid_z(valid)    = zscore(evid(valid));
                rt_z(valid)      = zscore(rt(valid));

            end
        end

        keep = rt >= 0 & ...
            ~isnan(evid_z)    & ...
            ~isnan(AvgWin_z)  & ...
            ~isnan(AvgLose_z) & ...
            ~isnan(conf)      & ...
            ~isnan(rt_z);

        %% ----------------------------------------------------
        % CONFIDENCE REGRESSION
        %% ----------------------------------------------------

        mdl_conf = fitglm( ...
            table( ...
            evid_z(keep), ...
            AvgWin_z(keep), ...
            AvgLose_z(keep), ...
            conf(keep)), ...
            'Var4 ~ Var1 + Var2 + Var3', ...
            'Distribution', 'binomial');

        betas_conf_time{tt} = [ ...
            betas_conf_time{tt}; ...
            mdl_conf.Coefficients.Estimate(2:4)'];

        %% ----------------------------------------------------
        % RT REGRESSION
        %% ----------------------------------------------------

        mdl_rt = fitlm( ...
            table( ...
            evid_z(keep), ...
            AvgWin_z(keep), ...
            AvgLose_z(keep), ...
            rt_z(keep)), ...
            'Var4 ~ Var1 + Var2 + Var3');

        betas_rt_time{tt} = [ ...
            betas_rt_time{tt}; ...
            mdl_rt.Coefficients.Estimate(2:4)'];

    end
end

%% ============================================================
% PLOT ? CONFIDENCE
%% ============================================================

figure('Color','w','Position',[100 100 900 350]);

cols = [ ...
    0 0 0;          % AvgWinDV
    0.00 0.45 0.90  % AvgLoseDV
];

Xw = nan(nTime, nSubjects);
Xl = nan(nTime, nSubjects);

for tt = 1:nTime

    B = betas_conf_time{tt};

    if ~isempty(B)
        Xw(tt,1:size(B,1)) = B(:,2);
        Xl(tt,1:size(B,1)) = B(:,3);
    end
end

sig_w = cluster_perm_1samp(Xw, nPerm, 0.05);
sig_l = cluster_perm_1samp(Xl, nPerm, 0.05);

mw = mean(Xw, 2, 'omitnan');
sw = std(Xw, [], 2, 'omitnan') ./ sqrt(sum(~isnan(Xw),2));

ml = mean(Xl, 2, 'omitnan');
sl = std(Xl, [], 2, 'omitnan') ./ sqrt(sum(~isnan(Xl),2));


sig_x=max([abs(mw);abs(ml)])*1.5;
hold on

shadedErrorBar( ...
    time_sec, mw, sw, ...
    'lineprops', {'Color', cols(1,:), 'LineWidth', 2});

shadedErrorBar( ...
    time_sec, ml, sl, ...
    'lineprops', {'Color', cols(2,:), 'LineWidth', 2});

scatter( ...
    time_sec(sig_w), ...
    repmat(sig_x+0.02, sum(sig_w), 1), ...
    30, cols(1,:), 'filled');

scatter( ...
    time_sec(sig_l), ...
    repmat(sig_x, sum(sig_l), 1), ...
    30, cols(2,:), 'filled');

yline(0,'k--');
xline(0,'k--');
xline(1,'k--')

xlabel('Time from stimulus onset (s)');
ylabel('\beta');
title(['Confidence ' set_name]);

set(gca,'FontSize',14);
box off;

%% ============================================================
% PLOT ? REACTION TIME
%% ============================================================

figure('Color','w','Position',[100 100 900 350]);

Xw = nan(nTime, nSubjects);
Xl = nan(nTime, nSubjects);

for tt = 1:nTime

    B = betas_rt_time{tt};

    if ~isempty(B)
        Xw(tt,1:size(B,1)) = B(:,2);
        Xl(tt,1:size(B,1)) = B(:,3);
    end
end

sig_w = cluster_perm_1samp(Xw, nPerm, 0.05);
sig_l = cluster_perm_1samp(Xl, nPerm, 0.05);

mw = mean(Xw, 2, 'omitnan');
sw = std(Xw, [], 2, 'omitnan') ./ sqrt(sum(~isnan(Xw),2));

ml = mean(Xl, 2, 'omitnan');
sl = std(Xl, [], 2, 'omitnan') ./ sqrt(sum(~isnan(Xl),2));
sig_x=max([abs(mw);abs(ml)])*1.5;

hold on

shadedErrorBar( ...
    time_sec, mw, sw, ...
    'lineprops', {'Color', cols(1,:), 'LineWidth', 2});

shadedErrorBar( ...
    time_sec, ml, sl, ...
    'lineprops', {'Color', cols(2,:), 'LineWidth', 2});

scatter( ...
    time_sec(sig_w), ...
    repmat(sig_x+0.02, sum(sig_w), 1), ...
    30, cols(1,:), 'filled');

scatter( ...
    time_sec(sig_l), ...
    repmat(sig_x, sum(sig_l), 1), ...
    30, cols(2,:), 'filled');

yline(0,'k--');
xline(0,'k--');
xline(1,'k--');

xlabel('Time from stimulus onset (s)');
ylabel('\beta');
title(['Reaction Time ' set_name]);

set(gca,'FontSize',14);
box off;