%% ============================================================
% PMd-projected Correlation Matrix
% Best 25 ROIs + PMd/M1
% Last sample window (last 100 ms of form stimulus offset)
%% ============================================================

clear;
close all;
clc;
rng(42);

%% ============================================================
% SETTINGS
%% ============================================================

subjects = { ...
    'S01','S04','S06','S07','S12','S15','S17','S18', ...
    'S11','S05','S19','S20','S22','S10','S02','S08', ...
    'S24','S21','S23','S25'};

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
% TIME SETTINGS
%% ============================================================



% Final decision window used for correlation matrix
tt_idx = 141:151;

%% ============================================================
% AREA DEFINITIONS
%% ============================================================

% Best 25 ROIs
set2 = { ...
    'AIP','FEF','5mv','23c','i6-8','24dd','6mp','PFt','IP2','8Av', ...
    'p24pr','7PC','23d','PF','6ma','55b','7AL','24dv','DVT','PFm', ...
    '31a','RSC','5L','V6','8C'};

somato = {'4','3b','1','2','3a'}; 

% PMd/M1 reference areas
refAreas = {'6a','6d','4'};

% Only first area set
area_sets = {
    set2(~ismember(set2, refAreas))
};

set_names = {'All significant'}; 

%% ============================================================
% STORAGE
%% ============================================================

All_corr_sbj = [];
All_corr_sbj_res = [];

%% ============================================================
% MAIN LOOP
%% ============================================================

for s = 1   % only first area set

    areas = area_sets{s};

    for si = 1:numel(subjects)

        subj = subjects{si};
        fprintf('Processing %s\n', subj);

        subj_data   = struct();
        valid_areas = {};
        valid_ref   = {};

        %% ----------------------------------------------------
        % LOAD DATA
        %% ----------------------------------------------------

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

            % MATLAB-safe field name
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

        if isempty(valid_areas)
            continue
        end

        %% ----------------------------------------------------
        % FIND COMMON SESSIONS
        %% ----------------------------------------------------

        all_names = [valid_areas valid_ref];

        sess_all = cellfun( ...
            @(a) unique(subj_data.(a).sess), ...
            all_names, ...
            'UniformOutput', false);

        common_sess = sess_all{1};

        for i = 2:numel(sess_all)
            common_sess = intersect(common_sess, sess_all{i});
        end

        if isempty(common_sess)
            continue
        end

        %% ----------------------------------------------------
        % RESTRICT TO COMMON SESSIONS
        %% ----------------------------------------------------

        for i = 1:numel(all_names)

            A = all_names{i};

            keep = ismember( ...
                subj_data.(A).sess, ...
                common_sess);

            fn = fieldnames(subj_data.(A));

            for ff = 1:numel(fn)
                subj_data.(A).(fn{ff}) = ...
                    subj_data.(A).(fn{ff})(keep,:);
            end
        end

        %% ----------------------------------------------------
        % SHARED VARIABLES
        %% ----------------------------------------------------

        ref = valid_areas{1};

        choice     = subj_data.(ref).choice;
        sess       = subj_data.(ref).sess;
        choiceSign = sign(choice);

        %% ----------------------------------------------------
        % COMPUTE CORRELATION MATRIX
        %% ----------------------------------------------------

        [Corr_matrix Corr_matrix_res] = compute_timepointwise_corr_matrix( ...
            subj_data, ...
            valid_areas, ...
            valid_ref, ...
            tt_idx, ...
            choiceSign, ...
            sess);

        All_corr_sbj(:,:,si) = Corr_matrix;
        All_corr_sbj_res(:,:,si) = Corr_matrix_res;


    end
end

%% ============================================================
% GROUP-AVERAGE CORRELATION MATRIX
%% ============================================================

corrMat = mean(All_corr_sbj, 3, 'omitnan');
corrMat_res = mean(All_corr_sbj_res, 3, 'omitnan');

% Labels for plotting
area_labels = { ...
    'AIP','FEF','5mv','23c','i6-8','24dd','6mp','PFt','IP2','8Av', ...
    'p24pr','7PC','23d','PF','6ma','55b','7AL','24dv','DVT','PFm', ...
    '31a','RSC','5L','V6','8C','PMd/M1'};

%% ============================================================
% PLOT CORRELATION MATRIX
%% ============================================================

figure('Color','w');

subplot(1,2,1)
imagesc(corrMat, [-1 1]);
axis square;

% White background
set(gca, 'Color', [1 1 1]);

% Blue = negative, White = zero, Red = positive
colormap(flipud(brewermap([], 'RdBu')));

c = colorbar;
c.Label.String = 'Correlation coefficient';

set(gca, ...
    'XTick', 1:length(area_labels), ...
    'XTickLabel', area_labels, ...
    'YTick', 1:length(area_labels), ...
    'YTickLabel', area_labels, ...
    'TickLength', [0 0], ...
    'FontSize', 10);

xtickangle(45);

title('Correlation Matrix raw DV');

box off;



subplot(1,2,2)
imagesc(corrMat_res, [-1 1]);
axis square;

% White background
set(gca, 'Color', [1 1 1]);

% Blue = negative, White = zero, Red = positive
colormap(flipud(brewermap([], 'RdBu')));

c = colorbar;
c.Label.String = 'Correlation coefficient';

set(gca, ...
    'XTick', 1:length(area_labels), ...
    'XTickLabel', area_labels, ...
    'YTick', 1:length(area_labels), ...
    'YTickLabel', area_labels, ...
    'TickLength', [0 0], ...
    'FontSize', 10);

xtickangle(45);

title('Correlation Matrix residual DV');

box off;