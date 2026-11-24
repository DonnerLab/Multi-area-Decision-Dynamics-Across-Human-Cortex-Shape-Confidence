%% =======================================================================
% Figure 2D (PMd) and Figure S4C (M1)
%
% Variance of the bilateral choice decoder output (neural DV) over time
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

% Bilateral decoder output of the 180 areas (here 6a, 6d, 4):
% <subject>_<session>_<area>_CHOICEdecoder_alltrials_all.mat
decoderFolder = fullfile(projectFolder, 'data', 'decoder_output', 'ROI_bilateral_all_areas');

% Path to functions
addpath(fullfile(projectFolder, 'functions'));


%% =======================================================================
% SETTINGS
%% =======================================================================

subjects = { ...
    'S01','S04','S06','S07','S12','S15','S17','S18','S11','S05', ...
    'S19','S20','S22','S10','S02','S08','S24','S21','S23','S25'};

nSubjects = numel(subjects);

% ROIs: name and areas whose decoder output is averaged
roiNames = {'PMd', 'M1'};
roiAreas = {{'6a','6d'}, {'4'}};

% Non-overlapping windows of 5 samples (50 ms); time of the first sample
% of each window (s from stimulus onset)
nWindows     = 40;
windowLength = 5;
timeAxis     = -0.5:0.05:1.45;


%% =======================================================================
% LOOP OVER ROIs
%% =======================================================================

for roi = 1:numel(roiNames)

    % windows x subjects
    variance_DV = nan(nWindows,nSubjects);

    for subject = 1:nSubjects

        %% ---------------------------------------------------------------
        % LOAD DECODER OUTPUT (all sessions)
        %% ---------------------------------------------------------------

        files = dir(fullfile(decoderFolder, ...
            [subjects{subject} '_*_' roiAreas{roi}{1} '_CHOICEdecoder_alltrials_all.mat']));

        % windows x sessions x choice x stimulus category
        variance_session = nan(nWindows,numel(files),2,2);

        for file = 1:numel(files)

            session = sscanf(files(file).name, [subjects{subject} '_%d_']);

            logOdds = 0;

            for area = 1:numel(roiAreas{roi})
                D = load(fullfile(decoderFolder, sprintf('%s_%d_%s_CHOICEdecoder_alltrials_all.mat', ...
                    subjects{subject}, session, roiAreas{roi}{area})));
                logOdds = logOdds + double(D.log_odds_all) / numel(roiAreas{roi});
            end

            choice   = double(D.Y_all(:));          % -1 = left, +1 = right
            stimulus = double(D.stimulus_all(:));   % -1 = left, +1 = right


            %% -----------------------------------------------------------
            % VARIANCE ACROSS TRIALS WITHIN EACH CHOICE x STIMULUS CATEGORY
            %% -----------------------------------------------------------

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
                        variance_session(window,file,choiceSide,stimulusSide) = ...
                            var(mean(logOdds(trials,windowSamples),2));
                    end

                end
            end

        end

        % Average across categories and sessions
        variance_DV(:,subject) = nanmean(nanmean(nanmean(variance_session,4),3),2);

    end


    %% ===================================================================
    % FIGURE (Figure 2D for PMd, Figure S4C for M1)
    %% ===================================================================

    figure('Color','w','Name',roiNames{roi});
    hold on

    shadedErrorBar( ...
        timeAxis, ...
        nanmean(variance_DV,2), ...
        nanstd(variance_DV,0,2) ./ sqrt(nSubjects), ...
        'lineProps', {'-','color','k'});

    xline(0, '-', 'Color', [0.5 0.5 0.5]);
    xline(1, '-', 'Color', [0.5 0.5 0.5]);
    for sample = 1:9
        xline(0.1*sample, '--', 'Color', [0.5 0.5 0.5]);
    end

    xticks([-0.5 0 1 1.5])
    xlabel('Time from stimulus onset (s)')
    ylabel('Variance (decoder output)')
    title(roiNames{roi})
    box off
    set(gca,'fontsize',15)

end
