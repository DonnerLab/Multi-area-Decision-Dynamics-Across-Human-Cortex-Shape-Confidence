%% =======================================================================
% Inter-hemispheric neural DV correlation with RT and Confidence
%
% Figure 2H amd S4E for Neural DV
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

% Stimulus window used for xcorr
timeWindow = 51:151;

%% =======================================================================
% Main loop
%% =======================================================================

for s = 1:nSubjects

    for ses = 1:nSessions

        %% ---------------------------------------------------------------
        % Initialize session variables
        %% ---------------------------------------------------------------

        All_RT_ses      = [];
        All_Conf_ses    = [];
        All_Perf_all    = [];
        All_Choices_ses = [];

        %% ---------------------------------------------------------------
        % File names
        %% ---------------------------------------------------------------

        leftFile = sprintf( ...
            '%s_choice_l_session%d_PMd_M1_decoder_output_100hz.mat', ...
            SubjectID_2{s}, ses);

        rightFile = sprintf( ...
            '%s_choice_r_session%d_PMd_M1_decoder_output_100hz.mat', ...
            SubjectID_2{s}, ses);

        if ~isfile(leftFile) || ~isfile(rightFile)
            continue
        end

        %% ---------------------------------------------------------------
        % Load decoder outputs
        %% ---------------------------------------------------------------

        load(leftFile)
        Perf_L = permute(Prob_Performance, [3,2,1]);

        load(rightFile)
        Perf_R = permute(Prob_Performance, [3,2,1]);

        % Convert choices: left = 0, right = 1
        All_Choices(All_Choices == -1) = 0;

        %% ---------------------------------------------------------------
        % Concatenate behavioral + neural variables
        %% ---------------------------------------------------------------

        All_RT_ses      = zscore(All_RT);
        All_Conf_ses    = All_Conf;
        All_Perf_all    = [Perf_L Perf_R];
        All_Choices_ses = All_Choices;

        % Split hemispheres again for clarity
        All_Perf_L = All_Perf_all(:,1:201);
        All_Perf_R = All_Perf_all(:,202:402);

        nTrials = size(All_Choices_ses,1);

        %% ---------------------------------------------------------------
        % Assign winning and losing hemisphere DVs
        %% ---------------------------------------------------------------

        DV_win  = [];
        DV_lose = [];

        for tr = 1:nTrials

            if All_Choices_ses(tr) == 0
                % left choice
                DV_win(tr,:)  = All_Perf_L(tr,:);
                DV_lose(tr,:) = All_Perf_R(tr,:);
            else
                % right choice
                DV_win(tr,:)  = All_Perf_R(tr,:);
                DV_lose(tr,:) = All_Perf_L(tr,:);
            end
        end

        %% ---------------------------------------------------------------
        % Trial-wise inter-hemispheric cross-correlation
        %% ---------------------------------------------------------------

        All_corr = [];

        for tr = 1:nTrials

            All_corr(:,tr) = xcorr( ...
                zscore(DV_win(tr,timeWindow))', ...
                zscore(DV_lose(tr,timeWindow))', ...
                lag, ...
                "normalized");
        end

        % Average across lags (-100 to 100 ms)
        All_corr_mean = mean(All_corr(1:21,:),1);

        %% ---------------------------------------------------------------
        % Regression with RT
        %% ---------------------------------------------------------------

        corr_RT = glmfit( ...
            zscore(All_corr_mean)', ...
            All_RT_ses);

        %% ---------------------------------------------------------------
        % Regression with confidence
        %% ---------------------------------------------------------------

        corr_conf = glmfit( ...
            zscore(All_corr_mean'), ...
            All_Conf_ses - 1, ...
            'binomial', ...
            'link', ...
            'logit');

        %% ---------------------------------------------------------------
        % Store beta coefficients
        %% ---------------------------------------------------------------

        Corr_RT(s,ses)   = corr_RT(2);
        Corr_Conf(s,ses) = corr_conf(2);
    end
end

%% =======================================================================
% Average across sessions
%% =======================================================================

Corr_Conf = mean(Corr_Conf,2);
Corr_RT   = mean(Corr_RT,2);

%% =======================================================================
% Plot
%% =======================================================================

close all
figure(1)

%% -----------------------------------------------------------------------
% Confidence
%% -----------------------------------------------------------------------

subplot(1,2,1)

bar(1, mean(Corr_Conf), ...
    'facecolor', 'k', ...
    'facealpha', 0.5)

hold on

for s = 1:nSubjects

    scatter( ...
        1 + (rand(1)*2 - 1)*0.05, ...
        Corr_Conf(s), ...
        100, ...
        'k', ...
        'filled', ...
        'markerfacealpha', 0.5)
end

pval_conf = permutationTest( ...
    Corr_Conf, ...
    zeros(nSubjects,1), ...
    10000);

title(num2str(pval_conf))

ylim([-0.15 0.15])
ylabel('weight for  confidence report')
xticks([])

%% -----------------------------------------------------------------------
% Reaction time
%% -----------------------------------------------------------------------

subplot(1,2,2)

bar(1, mean(Corr_RT), ...
    'facecolor', 'k', ...
    'facealpha', 0.5)

hold on

for s = 1:nSubjects

    scatter( ...
        1 + (rand(1)*2 - 1)*0.05, ...
        Corr_RT(s), ...
        100, ...
        'k', ...
        'filled', ...
        'markerfacealpha', 0.5)
end

pval_rt = permutationTest( ...
    Corr_RT, ...
    zeros(nSubjects,1), ...
    10000);

title(num2str(pval_rt))

ylim([-0.15 0.15])
ylabel('correlation with reaction time')
xticks([])

set(gca,'fontsize',15)
box off