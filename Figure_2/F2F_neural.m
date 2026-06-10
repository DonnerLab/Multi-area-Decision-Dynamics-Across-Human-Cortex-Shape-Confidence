%% =======================================================================
% Inter-hemispheric correlation of winning vs losing DV
% Figure 2F for neural DV
%% =======================================================================

clear
close all

%% =======================================================================
% Paths
%% =======================================================================

% Folder where this script is located
scriptFolder = fileparts(mfilename('fullpath'));

% Go up two folders
projectFolder = fileparts(fileparts(scriptFolder));

% Data folder
hemiFolder = fullfile(projectFolder, 'data', 'decoder_output', 'hemi');

cd(hemiFolder)

%% =======================================================================
% Subject information
%% =======================================================================

SubjectID_2 = { ...
    'S01','S02','S04','S05','S07','S08','S11','S12','S15','S17', ...
    'S18','S10','S06','S19','S20','S21','S22','S23','S24','S25'};

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
% Main analysis
%% =======================================================================

for s = 1:20

    allChoices = [];
    allSides   = [];
    DV_win_all = [];
    DV_lose_all = [];

    for ses = 1:6

        sessionID = sessions(s, ses);

        leftFile = [SubjectID_2{s}, ...
            '_choice_l_session', num2str(sessionID), ...
            '_PMd_M1_decoder_output_100hz.mat'];

        rightFile = [SubjectID_2{s}, ...
            '_choice_r_session', num2str(sessionID), ...
            '_PMd_M1_decoder_output_100hz.mat'];

        if ~isfile(leftFile) || ~isfile(rightFile)
            continue
        end

        %% Load left and right decoder outputs

        load(leftFile)
        Perf_L = Prob_Performance;

        load(rightFile)
        Perf_R = Prob_Performance;

        %% Assign winning and losing DVs

        DV_win  = [];
        DV_lose = [];

        for tr = 1:size(All_Choices, 1)

            if All_Choices(tr) == -1
                % left choice - right hemisphere wins
                DV_win(tr, :)  = permute(Perf_R(1, :, tr), [3 2 1]);
                DV_lose(tr, :) = permute(Perf_L(1, :, tr), [3 2 1]);
            else
                % right choice - left hemisphere wins
                DV_win(tr, :)  = permute(Perf_L(1, :, tr), [3 2 1]);
                DV_lose(tr, :) = permute(Perf_R(1, :, tr), [3 2 1]);
            end

        end

        %% Stimulus evidence and stimulus side

        Evidence = [ ...
            mean(All_samples(:, 1:10), 2), ...
            mean(All_samples(:,11:20), 2)];

        % 0 = left stronger, 1 = right stronger
        stimulusSide = double((Evidence(:,2) - Evidence(:,1)) > 0);

        % Remove equal-evidence trials
        id_equal = abs(Evidence(:,1) - Evidence(:,2)) < 1e-4;
        stimulusSide(id_equal) = 3;

        %% Convert choices: left = 0, right = 1

        All_Choices(All_Choices == -1) = 0;

        %% Concatenate across sessions

        allChoices   = [allChoices; All_Choices];
        allSides     = [allSides; stimulusSide];
        DV_win_all   = [DV_win_all; DV_win];
        DV_lose_all  = [DV_lose_all; DV_lose];

    end

    %% ===============================================================
    % Correlation between winning and losing DV
    % measured around 0.95 s after stimulus onset
    %% ===============================================================

    for stimSide = 1:2
        for ch = 1:2

            trialID = find( ...
                allChoices == (ch - 1) & ...
                allSides   == (stimSide - 1));

            if isempty(trialID)
                Hemi_corr(s, ch, stimSide) = NaN;
                continue
            end

            DV_win_mean = mean(DV_win_all(trialID, 141:151), 2);
            DV_lose_mean = mean(DV_lose_all(trialID, 141:151), 2);

            r = corr( ...
                zscore(DV_win_mean), ...
                zscore(DV_lose_mean));

            Hemi_corr(s, ch, stimSide) = r;

        end
    end
end

%% =======================================================================
% Plot
%% =======================================================================

close all
figure

All_Corr_hemi = squeeze(mean(mean(Hemi_corr, 3, 'omitnan'), 2, 'omitnan'));

bar(1, mean(All_Corr_hemi, 'omitnan'), ...
    'FaceColor', 'k', ...
    'FaceAlpha', 0.5)

hold on

for s = 1:20
    scatter( ...
        1 + (rand(1)*2 - 1)*0.05, ...
        All_Corr_hemi(s), ...
        100, ...
        'k', ...
        'filled', ...
        'MarkerFaceAlpha', 0.5)
end

xticks([])
ylabel('Correlation coefficient (neural DV)')
set(gca, 'FontSize', 15)
box off
axis square