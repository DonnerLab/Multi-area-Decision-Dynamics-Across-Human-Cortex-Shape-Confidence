%% Clean Generalization Matrix Analysis
% Plot cross-temporal generalization matrix

clear; close all; clc;

%% ================================================================
%% PATHS
%% ================================================================

%  path generation
scriptFolder  = fileparts(mfilename('fullpath'));
projectFolder = fileparts(fileparts(scriptFolder));

matrix_path = fullfile(projectFolder, 'data', 'generalization matrix');
cd(matrix_path)

%% ================================================================
%% SUBJECT IDS
%% ================================================================

subject_ids = { ...
    'S01','S04','S06','S07','S12','S15','S17','S18', ...
    'S11','S05','S19','S20','S22','S10','S02','S08', ...
    'S24','S21','S23','S25'};

n_subjects = 20;
n_sessions = 6;

%% ================================================================
%% LOAD GENERALIZATION MATRICES
%% ================================================================

all_subjects = [];

for s = 1:n_subjects

    session_data = [];
    matrix_size = [];

    for ses = 1:n_sessions

        file_name = sprintf('%s_CHOICE_decoder_session%d_PCA_PMd_M1_Generalization_matrix.mat', ...
            subject_ids{s}, ses);

        if isfile(file_name)
            load(file_name, 'Decoder_Performance')

            session_data(:,:,ses) = Decoder_Performance;
            matrix_size = size(Decoder_Performance);

        elseif ~isempty(matrix_size)
            % Preserve dimensions if one session is missing
            session_data(:,:,ses) = nan(matrix_size);
        end
    end

    % Average across sessions for this subject
    all_subjects(:,:,s) = nanmean(session_data, 3);
end

%% ================================================================
%% GROUP AVERAGE MATRIX
%% ================================================================

group_matrix = nanmean(all_subjects, 3);

%% ================================================================
%% PLOT
%% ================================================================

figure;

cmap = brewermap(180, '*RdBu');

imagesc(group_matrix)
colormap(cmap)
colorbar

xlim([1 801])
ylim([1 801])
axis square

caxis([0.42 0.58])

xticks([1 201 601 681])
yticks([1 201 601 681])

xticklabels({'-0.5','0','1','1.2'})
yticklabels({'-0.5','0','1','1.2'})

xlabel({'Testing time'; 'from stimulus onset (s)'})
ylabel({'Training time'; 'from stimulus onset (s)'})

set(gca, 'FontSize', 15)
