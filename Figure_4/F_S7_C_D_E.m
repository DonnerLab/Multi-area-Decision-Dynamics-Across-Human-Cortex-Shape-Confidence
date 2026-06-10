%% ============================================================
% Correlation analysis of Winning vs Losing DV
%
%% ============================================================
% Correlation analysis of Winning vs Losing DV
%
% Plot 1:
% Correlation between NON-ABSOLUTE winning/losing DV
% - before PMd/M1 removal  (raw)
% - after  PMd/M1 removal  (residual)
%
% Plot 2:
% Correlation between ABSOLUTE winning/losing DV
% - after PMd/M1 removal only
%
% Plot 3:
% Correlation between N winning areas and:
% - absolute winning DV
% - absolute losing DV
%
% Only first area set used:
% Best 25 ROIs
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

data_path = fullfile( ...
    projectFolder, ...
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
% AREA SET (BEST 25 ROIs)
%% ============================================================

areas = { ...
    'AIP','FEF','5mv','23c','i6_8','24dd','6mp','PFt','IP2', ...
    '8Av','p24pr','7PC','23d','PF','6ma','55b','7AL','24dv', ...
    'DVT','PFm','31a','RSC','5L','V6','8C'};

refAreas = {'6a','6d','4'};

%% ============================================================
% STORAGE
%% ============================================================

% Plot 1
Corr_raw   = cell(1, nTime);
Corr_resid = cell(1, nTime);

% Plot 2
Corr_abs   = cell(1, nTime);

% Plot 3
Corr_Nwin_Win  = cell(1, nTime);
Corr_Nwin_Lose = cell(1, nTime);

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
    % Load data
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
        SS = [];

        for f = 1:numel(files)

            D = load(fullfile(files(f).folder, files(f).name));

            X = D.log_odds_all(:,1:201);

            L  = [L; X];
            Y  = [Y; D.Y_all(:)];
            SS = [SS; f * ones(size(X,1),1)];

        end

        A = matlab.lang.makeValidName(area);

        subj_data.(A).logodds = L;
        subj_data.(A).choice  = Y;
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
    % Keep only common sessions across all areas
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
    % Shared variables
    %% --------------------------------------------------------

    ref = valid_areas{1};

    choice     = subj_data.(ref).choice;
    sess       = subj_data.(ref).sess;
    choiceSign = sign(choice);

    %% ========================================================
    % Sliding-window analysis
    %% ========================================================

    for tt = 1:nTime

        idx_t = win_starts(tt):(win_starts(tt) + win_length - 1);

        [AvgWinDV, AvgLoseDV, ...
         AvgWinDV_raw, AvgLoseDV_raw, ...
         AvgWinDV_nabs, AvgLoseDV_nabs, ...
         AvgWinDV_raw_nabs, AvgLoseDV_raw_nabs, ...
         N_win_areas] = ...
            compute_AvgWin_AvgLose_timepointwise_correlations( ...
                subj_data, ...
                valid_areas, ...
                valid_ref, ...
                idx_t, ...
                choiceSign, ...
                sess);

        %% ----------------------------------------------------
        % Session-wise z-scoring
        %% ----------------------------------------------------

        % Raw non-absolute
        Win_raw_z  = nan(size(AvgWinDV_raw_nabs));
        Lose_raw_z = nan(size(AvgLoseDV_raw_nabs));

        % Residual non-absolute
        Win_res_z  = nan(size(AvgWinDV_nabs));
        Lose_res_z = nan(size(AvgLoseDV_nabs));

        % Residual absolute
        Win_abs_z  = nan(size(AvgWinDV));
        Lose_abs_z = nan(size(AvgLoseDV));

        % Number of winning areas
        Nwin_z     = nan(size(N_win_areas));

        for ss = unique(sess)'

            idx = sess == ss;

            %% raw non-absolute
            valid = idx & ...
                ~isnan(AvgWinDV_raw_nabs) & ...
                ~isnan(AvgLoseDV_raw_nabs);

            if sum(valid) > 2
                Win_raw_z(valid)  = zscore(AvgWinDV_raw_nabs(valid));
                Lose_raw_z(valid) = zscore(AvgLoseDV_raw_nabs(valid));
            end

            %% residual non-absolute
            valid = idx & ...
                ~isnan(AvgWinDV_nabs) & ...
                ~isnan(AvgLoseDV_nabs);

            if sum(valid) > 2
                Win_res_z(valid)  = zscore(AvgWinDV_nabs(valid));
                Lose_res_z(valid) = zscore(AvgLoseDV_nabs(valid));
            end

            %% residual absolute
            valid = idx & ...
                ~isnan(AvgWinDV) & ...
                ~isnan(AvgLoseDV);

            if sum(valid) > 2
                Win_abs_z(valid)  = zscore(AvgWinDV(valid));
                Lose_abs_z(valid) = zscore(AvgLoseDV(valid));
            end

            %% number of winning areas
            valid = idx & ...
                ~isnan(N_win_areas);

            if sum(valid) > 2
                Nwin_z(valid) = zscore(N_win_areas(valid));
            end
        end

        %% ----------------------------------------------------
        % Plot 1:
        % Raw vs residual non-absolute correlation
        %% ----------------------------------------------------

        keep = ~isnan(Win_raw_z) & ~isnan(Lose_raw_z);

        if sum(keep) > 2
            Corr_raw{tt} = [ ...
                Corr_raw{tt}; ...
                corr(Win_raw_z(keep), Lose_raw_z(keep), ...
                'rows','complete')];
        end

        keep = ~isnan(Win_res_z) & ~isnan(Lose_res_z);

        if sum(keep) > 2
            Corr_resid{tt} = [ ...
                Corr_resid{tt}; ...
                corr(Win_res_z(keep), Lose_res_z(keep), ...
                'rows','complete')];
        end

        %% ----------------------------------------------------
        % Plot 2:
        % Absolute winning vs losing correlation
        %% ----------------------------------------------------

        keep = ~isnan(Win_abs_z) & ~isnan(Lose_abs_z);

        if sum(keep) > 2
            Corr_abs{tt} = [ ...
                Corr_abs{tt}; ...
                corr(Win_abs_z(keep), Lose_abs_z(keep), ...
                'rows','complete')];
        end

        %% ----------------------------------------------------
        % Plot 3:
        % N winning areas vs absolute winning / losing DV
        %% ----------------------------------------------------

        keep = ~isnan(Nwin_z) & ~isnan(Win_abs_z);

        if sum(keep) > 2
            Corr_Nwin_Win{tt} = [ ...
                Corr_Nwin_Win{tt}; ...
                corr(Nwin_z(keep), Win_abs_z(keep), ...
                'rows','complete')];
        end

        keep = ~isnan(Nwin_z) & ~isnan(Lose_abs_z);

        if sum(keep) > 2
            Corr_Nwin_Lose{tt} = [ ...
                Corr_Nwin_Lose{tt}; ...
                corr(Nwin_z(keep), Lose_abs_z(keep), ...
                'rows','complete')];
        end
    end
end

%% ============================================================
% Convert to matrices
%% ============================================================

X_raw        = nan(nTime, nSubjects);
X_resid      = nan(nTime, nSubjects);
X_abs        = nan(nTime, nSubjects);
X_Nwin_Win   = nan(nTime, nSubjects);
X_Nwin_Lose  = nan(nTime, nSubjects);

for tt = 1:nTime

    if ~isempty(Corr_raw{tt})
        X_raw(tt,1:numel(Corr_raw{tt})) = Corr_raw{tt};
    end

    if ~isempty(Corr_resid{tt})
        X_resid(tt,1:numel(Corr_resid{tt})) = Corr_resid{tt};
    end

    if ~isempty(Corr_abs{tt})
        X_abs(tt,1:numel(Corr_abs{tt})) = Corr_abs{tt};
    end

    if ~isempty(Corr_Nwin_Win{tt})
        X_Nwin_Win(tt,1:numel(Corr_Nwin_Win{tt})) = ...
            Corr_Nwin_Win{tt};
    end

    if ~isempty(Corr_Nwin_Lose{tt})
        X_Nwin_Lose(tt,1:numel(Corr_Nwin_Lose{tt})) = ...
            Corr_Nwin_Lose{tt};
    end
end

%% ============================================================
% Plot 1
%% ============================================================

figure('Color','w','Position',[100 100 900 350]);
hold on

shadedErrorBar( ...
    time_sec, ...
    mean(X_raw,2,'omitnan'), ...
    std(X_raw,[],2,'omitnan') ./ sqrt(sum(~isnan(X_raw),2)), ...
    'lineprops', {'k','LineWidth',2});

shadedErrorBar( ...
    time_sec, ...
    mean(X_resid,2,'omitnan'), ...
    std(X_resid,[],2,'omitnan') ./ sqrt(sum(~isnan(X_resid),2)), ...
    'lineprops', {'r','LineWidth',2});

yline(0,'k--');
xline(0,'k--');
xline(1,'k--');

xlabel('Time (s)');
ylabel('Correlation');
title('correlation between DVs of winning and losing areas');

legend({'Before removal of PMd/M1','After removal of PMd/M1'}, ...
    'Location','best');
ylim([-1 0])
set(gca,'FontSize',14);
box off;

%% ============================================================
% Plot 2
%% ============================================================

figure('Color','w','Position',[100 100 900 350]);
hold on

shadedErrorBar( ...
    time_sec, ...
    mean(X_abs,2,'omitnan'), ...
    std(X_abs,[],2,'omitnan') ./ sqrt(sum(~isnan(X_abs),2)), ...
    'lineprops', {'b','LineWidth',2});

yline(0,'k--');
xline(0,'k--');
xline(1,'k--');

xlabel('Time (s)');
ylabel('Correlation');
title('rectified winning vs rectified losing DV correlation');


ylim([-0.2 0.15])

set(gca,'FontSize',14);
box off;

%% ============================================================
% Plot 3
%% ============================================================

figure('Color','w','Position',[100 100 900 350]);
hold on

shadedErrorBar( ...
    time_sec, ...
    mean(X_Nwin_Win,2,'omitnan'), ...
    std(X_Nwin_Win,[],2,'omitnan') ./ sqrt(sum(~isnan(X_Nwin_Win),2)), ...
    'lineprops', {'k','LineWidth',2});

shadedErrorBar( ...
    time_sec, ...
    mean(X_Nwin_Lose,2,'omitnan'), ...
    std(X_Nwin_Lose,[],2,'omitnan') ./ sqrt(sum(~isnan(X_Nwin_Lose),2)), ...
    'lineprops', {'r','LineWidth',2});

yline(0,'k--');
xline(0,'k--');
xline(1,'k--');

xlabel('Time (s)');
ylabel('Correlation');
title('N winning areas vs absolute DV');

legend({ ...
    'N winning areas vs AvgWinDV', ...
    'N winning areas vs AvgLoseDV'}, ...
    'Location','best');

set(gca,'FontSize',14);
box off;