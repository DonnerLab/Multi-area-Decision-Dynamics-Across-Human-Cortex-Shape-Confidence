%% Variance and mean of decoder output 
clear;
close all;
clc;

%% =======================================================================
% Paths
%% =======================================================================

% Folder where this script is located
scriptFolder = fileparts(mfilename('fullpath'));

% Go up two folders
projectFolder = fileparts(fileparts(scriptFolder));

% Main folders
data_path       = fullfile(projectFolder, 'data', 'decoder_output', 'hemi');
functionsFolder = fullfile(projectFolder, 'functions');
giftiFolder     = fullfile(functionsFolder, 'gifti-1.6');

% Add required paths
addpath(data_path)
addpath(functionsFolder)
addpath(giftiFolder)

% Move to decoder folder
cd(data_path)

%% =======================================================================
% Subject information
%% =======================================================================

SubjectID_2 = { ...
    'S01','S02','S04','S05','S07','S08','S11','S12', ...
    'S15','S17','S18','S10','S06','S19','S20','S21', ...
    'S22','S23','S24','S25'};

Subjects = [ ...
    1,2,4,5,7,8,11,12,15,17,18,10,6,19,20,21,22,23,24,25];

Drugs = [ ...
    3 1 1 2 3 2;
    3 1 3 2 1 2;
    1 3 2 1 3 2;
    2 1 1 3 2 3;
    1 3 2 3 1 2;
    3 1 1 2 2 3;
    1 2 2 3 3 1;
    3 1 2 1 3 2;
    1 3 3 2 2 1;
    2 1 3 2 3 1;
    1 2 3 2 1 3;
    1 2 2 3 1 3;
    2 1 3 1 2 3;
    2 3 1 2 1 3;
    1 2 3 1 3 2;
    2 1 1 2 3 3;
    3 1 2 2 1 3;
    2 3 2 1 1 3;
    3 1 1 2 3 2;
    3 2 1 3 1 2];

Extra_sessions = [ ...
    1,2,4,5,7,0,11,12,15,17,18,10,0,0,0,0,0,0,0,0,0];

hemis = {'l','r'};
choices = [1, -1];

nSubjects = 20;
nSessions = 6;
time_axis = -0.5:0.05:1.45;

%% =======================================================================
% Main loop
%% =======================================================================

for s = 1:nSubjects

    session_list = 1:6;

    for ses = 1:nSessions

        file_left = [ ...
            SubjectID_2{s}, ...
            '_choice_l_session', ...
            num2str(session_list(ses)), ...
            '_PMd_M1_decoder_output_100hz.mat'];

        file_right = [ ...
            SubjectID_2{s}, ...
            '_choice_r_session', ...
            num2str(session_list(ses)), ...
            '_PMd_M1_decoder_output_100hz.mat'];

        if exist(file_left, 'file') ~= 2 || exist(file_right, 'file') ~= 2
            continue
        end

        %% Load left hemisphere
        load(file_left)
        Perf_L = Prob_Performance;

        %% Load right hemisphere
        load(file_right)
        Perf_R = Prob_Performance;

        %% Build contra / ipsi DV relative to choice

        DV_contra = [];
        DV_ipsi   = [];

        for tr = 1:size(All_Choices,1)

            if All_Choices(tr) == -1
                % left choice  contra = right hemi
                DV_contra(tr,:) = permute(Perf_R(1,:,tr), [3 2 1]);
                DV_ipsi(tr,:)   = permute(Perf_L(1,:,tr), [3 2 1]);

            else
                % right choice  contra = left hemi
                DV_contra(tr,:) = permute(Perf_L(1,:,tr), [3 2 1]);
                DV_ipsi(tr,:)   = permute(Perf_R(1,:,tr), [3 2 1]);
            end
        end

        %% Stimulus side
        Sides = sign( ...
            mean(All_samples(:,1:10),2) - ...
            mean(All_samples(:,11:20),2));

        % convert choices from {-1,1} to {0,1}
        All_Choices(All_Choices == -1) = 0;

        %% Loop across choice and stimulus side

        for ch = 1:2
            for sid = 1:2

                % trials for this choice + stimulus category
                id = find( ...
                    Sides == sid-1 & ...
                    All_Choices == ch-1);

                if isempty(id)
                    continue
                end

                %% Across time windows

                for time_p = 1:40

                    win = (1:5) + 5*(time_p-1);

                    %% Variance

                    var_contra = var( ...
                        mean(DV_contra(id, win), 2));

                    var_ipsi = var( ...
                        mean(DV_ipsi(id, win), 2));

                    %% Mean

                    mean_contra = mean( ...
                        mean(DV_contra(id, win), 2));

                    mean_ipsi = mean( ...
                        mean(DV_ipsi(id, win), 2));

                    %% Fano factor

                    ff_contra = var_contra ./ ...
                        mean(mean(DV_contra(:, win), 2));

                    ff_ipsi = var_ipsi ./ ...
                        mean(mean(DV_ipsi(:, win), 2));

                    %% Store temporary values

                    Var_contra(time_p)  = var_contra;
                    Var_ipsi(time_p)    = var_ipsi;

                    Mean_contra(time_p) = mean_contra;
                    Mean_ipsi(time_p)   = mean_ipsi;

                    FF_contra(time_p)   = ff_contra;
                    FF_ipsi(time_p)     = ff_ipsi;
                end

                %% Point process variance (PPV)

                PPV_contra = min(FF_contra);
                PPV_ipsi   = min(FF_ipsi);

                %% VarCE

                All_VarCE_contra(:,ses,s,ch,sid) = ...
                    Var_contra - Mean_contra .* PPV_contra;

                All_VarCE_ipsi(:,ses,s,ch,sid) = ...
                    Var_ipsi - Mean_ipsi .* PPV_ipsi;

                %% Raw variance

                All_Var_contra(:,ses,s,ch,sid) = Var_contra;
                All_Var_ipsi(:,ses,s,ch,sid)   = Var_ipsi;

                %% Mean DV

                All_Mean_contra(:,ses,s,ch,sid) = Mean_contra;
                All_Mean_ipsi(:,ses,s,ch,sid)   = Mean_ipsi;

            end
        end
    end
end

%% =======================================================================
% Correct trials only
%% =======================================================================

% Correct:
% choice left + stimulus left
% choice right + stimulus right

All_Variance_contra_corr = mean(cat(4, ...
    mean(mean(All_Var_contra(:,:,:,1,1),5),4), ...
    mean(mean(All_Var_contra(:,:,:,2,2),5),4)), 4);

All_Variance_ipsi_corr = mean(cat(4, ...
    mean(mean(All_Var_ipsi(:,:,:,1,1),5),4), ...
    mean(mean(All_Var_ipsi(:,:,:,2,2),5),4)), 4);

All_Mean_contra_corr = mean(cat(4, ...
    mean(mean(All_Mean_contra(:,:,:,1,1),5),4), ...
    mean(mean(All_Mean_contra(:,:,:,2,2),5),4)), 4);

All_Mean_ipsi_corr = mean(cat(4, ...
    mean(mean(All_Mean_ipsi(:,:,:,1,1),5),4), ...
    mean(mean(All_Mean_ipsi(:,:,:,2,2),5),4)), 4);

%% =======================================================================
% Plot variance
%% =======================================================================

figure(1);
clf;

subplot(1,2,1)

shadedErrorBar( ...
    time_axis, ...
    nanmean(nanmean(All_Variance_contra_corr,3),2), ...
    nanstd(nanmean(All_Variance_contra_corr,2),0,3) ...
    ./ sqrt(size(All_Variance_contra_corr,3)), ...
    'lineProps', {'-','color','k'});

hold on

shadedErrorBar( ...
    time_axis, ...
    nanmean(nanmean(All_Variance_ipsi_corr,3),2), ...
    nanstd(nanmean(All_Variance_ipsi_corr,2),0,3) ...
    ./ sqrt(size(All_Variance_ipsi_corr,3)), ...
    'lineProps', {'-','color','b'});

xticks([1 8 24 30])
xticklabels([-0.5 0 1 1.5])

xlabel('Time from stimulus onset')
ylabel('Variance (decoder output)')

set(gca,'fontsize',15)
box off

%% =======================================================================
% Plot mean neural DV
%% =======================================================================

subplot(1,2,2)

shadedErrorBar( ...
    time_axis, ...
    nanmean(nanmean(All_Mean_contra_corr,3),2), ...
    nanstd(nanmean(All_Mean_contra_corr,2),0,3) ...
    ./ sqrt(size(All_Mean_contra_corr,3)), ...
    'lineProps', {'-','color','k'});

hold on

shadedErrorBar( ...
    time_axis, ...
    nanmean(nanmean(All_Mean_ipsi_corr,3),2), ...
    nanstd(nanmean(All_Mean_ipsi_corr,2),0,3) ...
    ./ sqrt(size(All_Mean_ipsi_corr,3)), ...
    'lineProps', {'-','color','b'});

xticks([1 8 24 30])
xticklabels([-0.5 0 1 1.5])

xlabel('Time from stimulus onset')
ylabel('Neural DV (log odds)')

set(gca,'fontsize',15)
box off