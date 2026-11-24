%% =======================================================================
% Figure 2I (PMd) and Figure S6F (M1)
%
% Within-trial correlation between the winning and losing neural DV
% (single-hemisphere choice decoders), shown as cross-correlograms for
% two trial segments:
%
%   pre-stimulus interval : -0.5 - 0 s
%   early stimulus        :  0 - 0.5 s
%   and their difference (early - pre-stimulus, red)
%
% Each trial: linear detrending and z-scoring of the two DVs, normalized
% cross-correlation (lags +-100 ms)
% PMd = decoder output of areas 6a and 6d averaged; M1 = area 4
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

% Single-hemisphere decoder output
decoderFolder = fullfile(projectFolder, 'data', 'decoder_output', 'PMd_M1_hemi');

% Path to functions
addpath(fullfile(projectFolder, 'functions'));


%% =======================================================================
% SETTINGS
%% =======================================================================

subjects = { ...
    'S01','S02','S04','S05','S07','S08','S11','S12','S15','S17', ...
    'S18','S10','S06','S19','S20','S21','S22','S23','S24','S25'};

nSubjects     = numel(subjects);
nPermutations = 10000;    % permutations of the permutation tests

% ROIs: name and areas whose decoder output is averaged
roiNames = {'PMd', 'M1'};
roiAreas = {{'6a','6d'}, {'4'}};

% Trial segments (samples; 10 ms each, first sample at -0.5 s)
preStimulusSamples   = 1:51;     % -0.5 - 0 s
earlyStimulusSamples = 51:101;   %  0 - 0.5 s

% Cross-correlation lags: +-10 samples = +-100 ms
maxLag  = 10;
lagAxis = (-maxLag:maxLag) * 10;
nLags   = numel(lagAxis);


%% =======================================================================
% LOOP OVER ROIs
%% =======================================================================

for roi = 1:numel(roiNames)

    % subjects x lags
    crossCorrelation_pre   = nan(nSubjects,nLags);
    crossCorrelation_early = nan(nSubjects,nLags);

    for subject = 1:nSubjects

        files = dir(fullfile(decoderFolder, ...
            [subjects{subject} '_SES*_l_CHOICEdecoder_' roiAreas{roi}{1} '_all.mat']));

        % lags x sessions
        crossCorrelation_pre_session   = nan(nLags,numel(files));
        crossCorrelation_early_session = nan(nLags,numel(files));

        for file = 1:numel(files)

            %% -----------------------------------------------------------
            % LOAD THE DVs OF BOTH HEMISPHERES (this session)
            %% -----------------------------------------------------------

            session = sscanf(files(file).name, [subjects{subject} '_SES%d']);

            areaDV_right = [];   % left hemisphere
            areaDV_left  = [];   % right hemisphere

            for area = 1:numel(roiAreas{roi})

                L = load(fullfile(decoderFolder, sprintf('%s_SES%d_l_CHOICEdecoder_%s_all.mat', ...
                    subjects{subject}, session, roiAreas{roi}{area})));
                R = load(fullfile(decoderFolder, sprintf('%s_SES%d_r_CHOICEdecoder_%s_all.mat', ...
                    subjects{subject}, session, roiAreas{roi}{area})));

                areaDV_right(:,:,area) = double(L.log_odds_all);
                areaDV_left(:,:,area)  = double(R.log_odds_all);

            end

            DV_right = mean(areaDV_right,3);
            DV_left  = mean(areaDV_left,3);

            % Winning DV = DV of the chosen side (Y_all: +1 = right choice)
            rightChoice = double(L.Y_all(:)) == 1;

            winningDV = DV_left;
            losingDV  = DV_right;
            winningDV(rightChoice,:) = DV_right(rightChoice,:);
            losingDV(rightChoice,:)  = DV_left(rightChoice,:);


            %% -----------------------------------------------------------
            % SINGLE-TRIAL CROSS-CORRELATIONS
            %% -----------------------------------------------------------

            nTrials = size(winningDV,1);

            trialCorrelation_pre   = nan(nLags,nTrials);
            trialCorrelation_early = nan(nLags,nTrials);

            for trial = 1:nTrials

                trialCorrelation_pre(:,trial) = xcorr( ...
                    zscore(detrend(winningDV(trial,preStimulusSamples)))', ...
                    zscore(detrend(losingDV(trial,preStimulusSamples)))', ...
                    maxLag, 'normalized');

                trialCorrelation_early(:,trial) = xcorr( ...
                    zscore(detrend(winningDV(trial,earlyStimulusSamples)))', ...
                    zscore(detrend(losingDV(trial,earlyStimulusSamples)))', ...
                    maxLag, 'normalized');

            end

            crossCorrelation_pre_session(:,file)   = mean(trialCorrelation_pre,2);
            crossCorrelation_early_session(:,file) = mean(trialCorrelation_early,2);

        end

        % Average across sessions
        crossCorrelation_pre(subject,:)   = nanmean(crossCorrelation_pre_session,2);
        crossCorrelation_early(subject,:) = nanmean(crossCorrelation_early_session,2);

    end

    crossCorrelation_difference = crossCorrelation_early - crossCorrelation_pre;


    %% ===================================================================
    % STATISTICS
     %% ===================================================================

    allCorrelograms = {crossCorrelation_pre, crossCorrelation_early, crossCorrelation_difference};

    significantLags = nan(3,nLags);

    for correlogram = 1:3

        pValues = nan(1,nLags);

        for lag = 1:nLags
            pValues(lag) = permutationTest( ...
                allCorrelograms{correlogram}(:,lag), ...
                zeros(nSubjects,1), ...
                nPermutations);
        end

        significantLags(correlogram,:) = double(fdr_bh(pValues,0.05));

    end

    significantLags(significantLags == 0) = nan;


    %% ===================================================================
    % FIGURE (Figure 2I for PMd, Figure S6F for M1)
    %
    % black: pre-stimulus, blue: early stimulus, red: early - pre
    %% ===================================================================

    correlogramColors = [0 0 0; 0 0 1; 0.8 0.2 0.2];

    figure('Color','w','Name',roiNames{roi});
    hold on

    for correlogram = 1:3

        shadedErrorBar( ...
            lagAxis, ...
            nanmean(allCorrelograms{correlogram},1), ...
            nanstd(allCorrelograms{correlogram},0,1) ./ sqrt(nSubjects), ...
            'lineProps', {'-','color',correlogramColors(correlogram,:)});

        % Significant lags (one row of marks per correlogram)
        scatter(lagAxis, significantLags(correlogram,:) .* (-0.17 - 0.03*correlogram), ...
            50, correlogramColors(correlogram,:), 'filled', 'square');

    end

    ylim([-0.29 0.01])
    xlabel('Time lag (ms)')
    ylabel('Correlation coefficient')
    title([roiNames{roi} ': pre-stimulus (black), early stimulus (blue), early - pre (red)'], ...
        'FontWeight', 'normal', 'FontSize', 10)
    box off
    set(gca,'fontsize',15)

end
