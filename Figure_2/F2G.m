%% =======================================================================
% Figure 2G
%
% Time course of the trial-averaged winning and losing neural DV
% (single-hemisphere choice decoders, PMd)
%
% PMd = decoder output of areas 6a and 6d averaged
%% =======================================================================

clear all
close all


%% =======================================================================
% PATHS
%% =======================================================================

% Folder where this script is located
scriptFolder = fileparts(mfilename('fullpath'));
if isempty(scriptFolder), scriptFolder = pwd; end

% Go up one folder (repository root)
projectFolder = fileparts(scriptFolder);

decoderFolder = fullfile(projectFolder, 'data', 'decoder_output', 'PMd_M1_hemi');

% Path to functions
addpath(fullfile(projectFolder, 'functions'));


%% =======================================================================
% SETTINGS
%% =======================================================================

subjects = { ...
    'S01','S02','S04','S05','S07','S08','S11','S12','S15','S17', ...
    'S18','S10','S06','S19','S20','S21','S22','S23','S24','S25'};

nSubjects = numel(subjects);

% Areas of PMd (decoder output averaged)
pmdAreas = {'6a','6d'};

% Non-overlapping windows of 5 samples (50 ms); time of the first sample
% of each window (s from stimulus onset)
nWindows     = 40;
windowLength = 5;
timeAxis     = -0.5:0.05:1.45;

% windows x subjects
mean_winningDV = nan(nWindows,nSubjects);
mean_losingDV  = nan(nWindows,nSubjects);


%% =======================================================================
% MAIN LOOP ACROSS SUBJECTS
%% =======================================================================

for subject = 1:nSubjects

    files = dir(fullfile(decoderFolder, ...
        [subjects{subject} '_SES*_l_CHOICEdecoder_' pmdAreas{1} '_all.mat']));

    % windows x sessions x choice x stimulus category
    mean_winning_session = nan(nWindows,numel(files),2,2);
    mean_losing_session  = nan(nWindows,numel(files),2,2);

    for file = 1:numel(files)

        %% ---------------------------------------------------------------
        % LOAD DECODER OUTPUT OF BOTH HEMISPHERES (this session)
        %
        % Labels from the left-hemisphere files:
        % Y_all: -1 = left choice, +1 = right choice
        % stimulus_all: -1 = left > right, +1 = right > left
        %% ---------------------------------------------------------------

        session = sscanf(files(file).name, [subjects{subject} '_SES%d']);

        DV_right = 0;
        DV_left  = 0;

        for area = 1:numel(pmdAreas)

            L = load(fullfile(decoderFolder, sprintf('%s_SES%d_l_CHOICEdecoder_%s_all.mat', ...
                subjects{subject}, session, pmdAreas{area})));
            R = load(fullfile(decoderFolder, sprintf('%s_SES%d_r_CHOICEdecoder_%s_all.mat', ...
                subjects{subject}, session, pmdAreas{area})));

            DV_right = DV_right + double(L.log_odds_all) / numel(pmdAreas);
            DV_left  = DV_left  + double(R.log_odds_all) / numel(pmdAreas);

        end

        choice   = double(L.Y_all(:));
        stimulus = double(L.stimulus_all(:));

        % Winning DV = DV of the chosen side, losing DV = DV of the other side
        winningDV = DV_left;
        losingDV  = DV_right;
        winningDV(choice == 1,:) = DV_right(choice == 1,:);
        losingDV(choice == 1,:)  = DV_left(choice == 1,:);


        %% ---------------------------------------------------------------
        % MEAN ACROSS TRIALS WITHIN EACH CHOICE x STIMULUS CATEGORY
        %% ---------------------------------------------------------------

        categoryValues = [-1 1];

        for choiceSide = 1:2
            for stimulusSide = 1:2

                trials = choice == categoryValues(choiceSide) & ...
                         stimulus == categoryValues(stimulusSide);

                if sum(trials) < 2
                    continue
                end

                for window = 1:nWindows
                    windowSamples = (1:windowLength) + windowLength*(window-1);
                    mean_winning_session(window,file,choiceSide,stimulusSide) = ...
                        mean(mean(winningDV(trials,windowSamples),2));
                    mean_losing_session(window,file,choiceSide,stimulusSide) = ...
                        mean(mean(losingDV(trials,windowSamples),2));
                end

            end
        end

    end

    % Average across categories and sessions
    mean_winningDV(:,subject) = nanmean(nanmean(nanmean(mean_winning_session,4),3),2);
    mean_losingDV(:,subject)  = nanmean(nanmean(nanmean(mean_losing_session,4),3),2);

end


%% =======================================================================
% FIGURE (Figure 2G)
%% =======================================================================

figure('Color','w');
hold on

shadedErrorBar( ...
    timeAxis, ...
    nanmean(mean_winningDV,2), ...
    nanstd(mean_winningDV,0,2) ./ sqrt(nSubjects), ...
    'lineProps', {'-','color','k'});

shadedErrorBar( ...
    timeAxis, ...
    nanmean(mean_losingDV,2), ...
    nanstd(mean_losingDV,0,2) ./ sqrt(nSubjects), ...
    'lineProps', {'-','color','b'});

yline(0,'--','color',[0.5 0.5 0.5]);

xticks([-0.5 0 1 1.5])
xlabel('Time from stimulus onset (s)')
ylabel('Neural DV (log odds)')
title('Winning DV (black), losing DV (blue)')
box off
set(gca,'fontsize',15)
