%% ============================================================
% PMd/M1 correlation with area decision variables (DV)
%
% 
% Figure S7A:
% Raw correlation before PMd/M1 projection removal
%   corr(PMd, DV_raw)
%
% Figure S7B:
% Residual correlation after PMd/M1 projection removal
%   corr(PMd, DV_residual)
%
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
% AREA DEFINITIONS
%% ============================================================

set2 = { ...
    'AIP','FEF','5mv','23c','i6-8','24dd','6mp','PFt','IP2', ...
    '8Av','p24pr','7PC','23d','PF','6ma','55b','7AL','24dv', ...
    'DVT','PFm','31a','RSC','5L','V6','8C'};

refAreas = {'6a','6d','4'};

% Use only first area set
areas = set2;

%% ============================================================
% STORAGE
%% ============================================================

% Raw correlation: corr(PMd, DV_raw)
Corr_raw = cell(1, nTime);

% Residual correlation: corr(PMd, DV_residual)
Corr_resid = cell(1, nTime);

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
    % LOAD DATA
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
    % KEEP ONLY COMMON SESSIONS
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
    % SHARED VARIABLES
    %% --------------------------------------------------------

    ref = valid_areas{1};

    choice     = subj_data.(ref).choice;
    sess       = subj_data.(ref).sess;
    choiceSign = sign(choice);

    %% ========================================================
    % SLIDING WINDOW
    %% ========================================================

    for tt = 1:nTime

        idx_t = win_starts(tt):(win_starts(tt) + win_length - 1);

        % correlation for this window
        [All_corr_raw, All_corr_resid] = ...
            compute_corr_DV_PMd_timepointwise( ...
                subj_data, ...
                valid_areas, ...
                valid_ref, ...
                idx_t, ...
                choiceSign, ...
                sess);

        Corr_raw{tt} = [ ...
            Corr_raw{tt}; ...
            All_corr_raw'];

        Corr_resid{tt} = [ ...
            Corr_resid{tt}; ...
            All_corr_resid'];

    end
end

%% ============================================================
% PLOT

%% ============================================================

figure( ...
    'Color', 'w', ...
    'Position', [100 100 1300 450]);

cols = brewermap(25, '*RdBu');

for ai = 1:25

    X_raw   = nan(nTime, nSubjects);
    X_resid = nan(nTime, nSubjects);

    for tt = 1:nTime

        B = [ ...
            Corr_raw{tt}(:, ai), ...
            Corr_resid{tt}(:, ai)];

        if ~isempty(B)

            X_raw(tt, 1:size(B,1)) = B(:,1);
            X_resid(tt, 1:size(B,1)) = B(:,2);

        end
    end

    %% --------------------------------------------------------
    % Mean + SEM
    %% --------------------------------------------------------

    mean_raw = mean(X_raw, 2, 'omitnan');
    sem_raw  = std(X_raw, [], 2, 'omitnan') ...
        ./ sqrt(sum(~isnan(X_raw), 2));

    mean_resid = mean(X_resid, 2, 'omitnan');
    sem_resid  = std(X_resid, [], 2, 'omitnan') ...
        ./ sqrt(sum(~isnan(X_resid), 2));

    %% --------------------------------------------------------
    % Plot raw
    %% --------------------------------------------------------

    subplot(1,2,1);
    hold on;

    shadedErrorBar( ...
        time_sec, ...
        mean_raw, ...
        sem_raw, ...
        'lineprops', { ...
            'Color', cols(ai,:), ...
            'LineWidth', 2});

    ylim([-0.05 1]);

    %% --------------------------------------------------------
    % Plot residual
    %% --------------------------------------------------------

    subplot(1,2,2);
    hold on;

    shadedErrorBar( ...
        time_sec, ...
        mean_resid, ...
        sem_resid, ...
        'lineprops', { ...
            'Color', cols(ai,:), ...
            'LineWidth', 2});

    ylim([-0.05 1]);

end

%% ============================================================
% FINAL FIGURE FORMATTING
%% ============================================================

subplot(1,2,1);

yline(0, 'k--');
xline(0, 'k--');
xline(0.5, 'k--');

xlabel('Time (s)');
ylabel('Correlation');

title('Raw correlation: corr(PMd, DV)');

set(gca, 'FontSize', 14);
box off;

subplot(1,2,2);

yline(0, 'k--');
xline(0, 'k--');
xline(0.5, 'k--');

xlabel('Time (s)');
ylabel('Correlation');

title('After PMd/M1 projection removal');

set(gca, 'FontSize', 14);
box off;