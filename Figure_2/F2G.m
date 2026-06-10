%% =======================================================================
% Inter-hemispheric cross-correlation of decoder output
%
% Figure 2G
%% =======================================================================

clear all
close all

%% =======================================================================
% Paths
%% =======================================================================

% Folder where this script is located
scriptFolder = fileparts(mfilename('fullpath'));

% Go up two folders
projectFolder = fileparts(fileparts(scriptFolder));

% Decoder output folder
hemiFolder = fullfile(projectFolder, 'data', 'decoder_output', 'hemi');

% Move to decoder folder
cd(hemiFolder)

%% =======================================================================
% Subject information
%% =======================================================================

SubjectID_2 = { ...
    'S01','S02','S04','S05','S07','S08','S11','S12','S15','S17', ...
    'S18','S10','S06','S19','S20','S21','S22','S23','S24','S25'};

SubjectID = [ ...
     1,4,6,7,12,15,17,18,11,05, ...
    19,20,22,10,02,08,24,21,23,25];

sessions = [ ...
    11,2,3,4,5,6;
    11,2,3,4,5,6;
     1,2,3,4,5,6;
    11,2,3,4,5,6;
    11,2,3,4,5,6;
    11,2,3,4,5,6;
    11,2,3,4,5,6;
    11,2,3,4,5,6;
    11,2,3,4,5,6;
    11,2,3,4,5,6;
     1,2,3,4,5,6;
     1,2,3,4,5,6;
     1,2,3,4,5,6;
    11,2,3,4,5,6;
    11,2,3,4,5,6;
     1,2,3,4,5,6;
     1,2,3,4,5,6;
     1,2,3,4,5,6;
     1,2,3,4,5,6;
     1,2,3,4,5,6];

%% =======================================================================
% Analysis parameters
%% =======================================================================

nSubjects = 20;
nSessions = 6;
lag       = 10;

% Time windows
idx_baseline = 1:51;
idx_early    = 51:101;
idx_late     = 101:151;
idx_post     = 151:201;

% Xcorr lags
lagAxis = -100:10:100;

%% =======================================================================
% Main loop
%% =======================================================================

for s = 1:nSubjects

    for ses = 1:nSessions

        sessionNumber = ses;

        leftFile = sprintf( ...
            '%s_choice_l_session%d_PMd_M1_decoder_output_100hz.mat', ...
            SubjectID_2{s}, sessionNumber);

        rightFile = sprintf( ...
            '%s_choice_r_session%d_PMd_M1_decoder_output_100hz.mat', ...
            SubjectID_2{s}, sessionNumber);

        if ~isfile(leftFile) || ~isfile(rightFile)
            continue
        end

        %% ---------------------------------------------------------------
        % Load decoder output
        %% ---------------------------------------------------------------

        load(leftFile)
        Perf_L = permute(Prob_Performance, [3,2,1]);

        load(rightFile)
        Perf_R = permute(Prob_Performance, [3,2,1]);

        % Convert choices: left = 0, right = 1
        All_Choices(All_Choices == -1) = 0;

        nTrials = size(All_Choices, 1);

        %% ---------------------------------------------------------------
        % Assign winning and losing hemisphere DVs
        %% ---------------------------------------------------------------

        DV_win  = [];
        DV_lose = [];

        for tr = 1:nTrials

            if All_Choices(tr) == 0
                % left choice - right hemi wins
                DV_win(tr,:)  = Perf_R(tr,:);
                DV_lose(tr,:) = Perf_L(tr,:);
            else
                % right choice - left hemi wins
                DV_win(tr,:)  = Perf_L(tr,:);
                DV_lose(tr,:) = Perf_R(tr,:);
            end
        end

        %% ---------------------------------------------------------------
        % Trial-wise cross-correlations
        %% ---------------------------------------------------------------

        All_corr_base  = [];
        All_corr_early = [];
        All_corr_late  = [];
        All_corr_post  = [];

        for tr = 1:nTrials

            % Baseline
            All_corr_base(:,tr) = xcorr( ...
                zscore(DV_win(tr,idx_baseline))', ...
                zscore(DV_lose(tr,idx_baseline))', ...
                lag, "normalized");

            % Early
            All_corr_early(:,tr) = xcorr( ...
                zscore(DV_win(tr,idx_early))', ...
                zscore(DV_lose(tr,idx_early))', ...
                lag, "normalized");

            % Late
            All_corr_late(:,tr) = xcorr( ...
                zscore(DV_win(tr,idx_late))', ...
                zscore(DV_lose(tr,idx_late))', ...
                lag, "normalized");

            % Post
            All_corr_post(:,tr) = xcorr( ...
                zscore(DV_win(tr,idx_post))', ...
                zscore(DV_lose(tr,idx_post))', ...
                lag, "normalized");
        end

        %% ---------------------------------------------------------------
        % Store average correlograms
        %% ---------------------------------------------------------------

        All_correlograms(:,:,s,1,ses) = [ ...
            mean(All_corr_base,2), ...
            mean(All_corr_early,2), ...
            mean(All_corr_late,2), ...
            mean(All_corr_post,2)];
    end
end

%% =======================================================================
% Plot baseline-corrected cross-correlograms
%% =======================================================================

colors = winter(6);
colors(1,:) = [0 0 0];

figure
hold on

for epoch = 2:4

    % Baseline correction
    All_corrgram = reshape( ...
        mean( ...
            All_correlograms(:,epoch,:,:,1) - ...
            All_correlograms(:,1,:,:,1), 4), ...
        21, 20);

    All_corrgram = All_corrgram';

    %% ---------------------------------------------------------------
    % Plot shaded error bar
    %% ---------------------------------------------------------------

    shadedErrorBar( ...
        lagAxis, ...
        nanmean(All_corrgram,1), ...
        std(All_corrgram)./sqrt(20), ...
        'lineProps', {'-','color',colors(epoch,:)});

    %% ---------------------------------------------------------------
    % Permutation test per lag
    %% ---------------------------------------------------------------

    P = nan(1,21);

    for ii = 1:21
        P(ii) = permutationTest( ...
            All_corrgram(:,ii), ...
            zeros(20,1), ...
            10000);
    end

    % FDR correction
    All_P(epoch,:) = double(fdr_bh(P,0.05));
    All_P(All_P == 0) = -10;

    %% ---------------------------------------------------------------
    % Significance markers
    %% ---------------------------------------------------------------

    scatter( ...
        lagAxis, ...
        All_P(epoch,:) .* (-0.13 - (0.01*(epoch-1))), ...
        50, ...
        colors(epoch,:), ...
        'filled', ...
        'square');
end

ylim([-0.17 0.01])

xlabel('Lag (ms)')
ylabel('Cross-correlation (baseline corrected)')
set(gca,'fontsize',15)
box off