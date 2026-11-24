%% =======================================================================
% Figure 3B-C, Figure S9B, C, F and Figure S10
%
% Regression of single-trial confidence (logistic) and RT (linear) on
% evidence strength and on the winning and losing neural DV
% (single-hemisphere choice decoders, winlose_regression.m)
%
%   left hemisphere decoder : DV for the right choice
%   right hemisphere decoder: DV for the left choice
%


%   PMd, confidence: Figure 3B (left), 3C (middle); right panel: Figure 3D
%                    (F3D_FS8B.m, with the model DVs)
%   PMd, RT        : Figure S9B (left), S9C (middle), S9F (right)
%   M1,  confidence: Figure S10A (middle), S10B (right)
%   M1,  RT        : Figure S10C (middle), S10D (right)
%
% PMd = decoder trained on areas 6a and 6d together ('6A6D'); M1 = area 4
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

% Single-hemisphere decoder output with confidence and RT of each trial:
decoderFolder = fullfile(projectFolder, 'data', 'decoder_output', 'PMd_M1_hemi_behavior');

% Path to functions
addpath(fullfile(projectFolder, 'functions'));


%% =======================================================================
% SETTINGS
%% =======================================================================

subjects = { ...
    'S01','S02','S04','S05','S06','S07','S08','S10','S11','S12', ...
    'S15','S17','S18','S19','S20','S21','S22','S23','S24','S25'};

nSubjects     = numel(subjects);
nPermutations = 10000;    % permutations of the permutation tests

% ROIs: name and decoder area
roiNames = {'PMd', 'M1'};
roiAreas = {'6A6D', '4'};

% Decoder time base: 201 samples (10 ms), first sample at -0.5 s
nTimePoints = 201;

% Sliding windows of 35 samples (350 ms), steps of 10 samples (100 ms)
windowLength = 35;
windowStarts = 1:10:(nTimePoints - windowLength + 1);
nWindows     = numel(windowStarts);
windowTime   = ((windowStarts + windowLength/2 - 1) * 10 - 500) / 1000;   % window centers (s)

% End of stimulus viewing
offsetWindow = 134:168;

behaviorNames = {'Confidence', 'RT'};



%% =======================================================================
% LOOP OVER ROIs
%% =======================================================================

for roi = 1:numel(roiNames)

    % subjects x windows x [evidence, winning DV, losing DV] x behavior
    betas_time   = nan(nSubjects,nWindows,3,2);

    % subjects x [evidence, winning DV, losing DV] x behavior
    betas_offset = nan(nSubjects,3,2);


    %% ===================================================================
    % MAIN LOOP ACROSS SUBJECTS
    %% ===================================================================

    for subject = 1:nSubjects

        %% ---------------------------------------------------------------
        % LOAD THE DVs OF BOTH HEMISPHERES AND THE BEHAVIOR (all sessions)
        %% ---------------------------------------------------------------

        files = dir(fullfile(decoderFolder, ...
            [subjects{subject} '_SES*_l_CHOICEdecoder_' roiAreas{roi} '_all.mat']));

        DV_R       = [];   % left hemisphere, DV for the right choice
        DV_L       = [];   % right hemisphere, DV for the left choice
        isRight    = [];   % choice (true = right)
        confidence = [];   % 1 = high, 0 = low
        RT         = [];
        evidence   = [];   % evidence strength
        session    = [];

        for file = 1:numel(files)

            sessionNumber = sscanf(files(file).name, [subjects{subject} '_SES%d']);

            L = load(fullfile(decoderFolder, sprintf('%s_SES%d_l_CHOICEdecoder_%s_all.mat', ...
                subjects{subject}, sessionNumber, roiAreas{roi})));
            R = load(fullfile(decoderFolder, sprintf('%s_SES%d_r_CHOICEdecoder_%s_all.mat', ...
                subjects{subject}, sessionNumber, roiAreas{roi})));

            DV_R       = [DV_R;       double(L.log_odds_all(:,1:nTimePoints))];
            DV_L       = [DV_L;       double(R.log_odds_all(:,1:nTimePoints))];
            isRight    = [isRight;    double(L.Y_all(:)) == 1];
            confidence = [confidence; double(L.confidence_all(:)) - 1];
            RT         = [RT;         double(L.rts_all(:))];
            evidence   = [evidence;   double(L.evidence_strength(:))];
            session    = [session;    sessionNumber * ones(size(L.log_odds_all,1),1)];

        end

        % Trials with a response
        keep = RT > 0;

        behavior = {confidence, RT};
        analysis = {'confidence', 'RT'};


        %% ---------------------------------------------------------------
        % REGRESSIONS: SLIDING WINDOWS AND END OF STIMULUS VIEWING
        %% ---------------------------------------------------------------

        for behaviorType = 1:2

            for window = 1:nWindows
                betas_time(subject,window,:,behaviorType) = winlose_regression( ...
                    DV_R, DV_L, isRight, behavior{behaviorType}, evidence, session, keep, ...
                    windowStarts(window):(windowStarts(window) + windowLength - 1), ...
                    analysis{behaviorType}, 0);
            end

            betas_offset(subject,:,behaviorType) = winlose_regression( ...
                DV_R, DV_L, isRight, behavior{behaviorType}, evidence, session, keep, ...
                offsetWindow, analysis{behaviorType}, 0);

        end

    end


    %% ===================================================================
    % FIGURES: ONE FOR EACH BEHAVIOR
    %% ===================================================================

    for behaviorType = 1:2

        figure('Color','w','Position',[100 100 1400 330], ...
            'Name',[roiNames{roi} ' - ' behaviorNames{behaviorType}]);

        % ---------------------------------------------------------------
        % Left: weight of evidence strength at the end of stimulus viewing
        % ---------------------------------------------------------------

        evidenceWeight = betas_offset(:,1,behaviorType);

        pValue_evidence = permutationTest(evidenceWeight, zeros(nSubjects,1), nPermutations);

        subplot(1,5,1)
        hold on

        bar(1, nanmean(evidenceWeight), 0.6, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'k');
        scatter(1 + (rand(nSubjects,1) - 0.5)*0.25, evidenceWeight, 70, 'k', 'filled', ...
            'MarkerFaceAlpha', 0.5);

        yline(0,'k-');
        xlim([0.4 1.6])
        xticks(1)
        xticklabels({'Evidence strength'})
        ylabel('Weight (a.u.)')
        title(sprintf('p = %.2g', pValue_evidence), 'FontWeight', 'normal')
        box off
        set(gca,'FontSize',14)

        % ---------------------------------------------------------------
        % Middle: weights of the winning (black) and losing (blue) DV
        % over time, with significant clusters
        % ---------------------------------------------------------------

        subplot(1,5,2:4)
        hold on

        dvColors = [0 0 0; 0.15 0.35 0.75];

        for dv = 1:2

            % windows x subjects (subjects with all windows fitted)
            weights = squeeze(betas_time(:,:,dv+1,behaviorType))';
            weights = weights(:, all(isfinite(weights),1));

            shadedErrorBar(windowTime, mean(weights,2), std(weights,0,2) ./ sqrt(size(weights,2)), ...
                'lineProps', {'-','Color',dvColors(dv,:),'LineWidth',2});

            % Cluster-based permutation test against zero
            significantWindows = double(cluster_perm_1samp(weights, nPermutations, 0.05));
            significantWindows(significantWindows == 0) = nan;

            meanWeights  = nanmean(betas_time(:,:,2:3,behaviorType),1);
            markerHeight = max(meanWeights(:)) * (1.3 - 0.15*dv);

            plot(windowTime, significantWindows * markerHeight, '-', ...
                'Color', dvColors(dv,:), 'LineWidth', 6);

        end

        yline(0,'k-');
        xline(0,'k--');
        xline(1,'k--');
        xlabel('Time from stimulus onset (s)')
        ylabel('Weight (a.u.)')
        title([roiNames{roi} ' neural DVs: winning (black), losing (blue)'], 'FontWeight', 'normal')
        box off
        set(gca,'FontSize',14)

        % ---------------------------------------------------------------
        % Right: |weight of winning DV| - |weight of losing DV| at the end
        % of stimulus viewing
        % ---------------------------------------------------------------

        weightDifference = abs(betas_offset(:,2,behaviorType)) - abs(betas_offset(:,3,behaviorType));

        pValue_difference = permutationTest(weightDifference, zeros(nSubjects,1), nPermutations);

        subplot(1,5,5)
        hold on

        bar(1, nanmean(weightDifference), 0.6, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'k');
        scatter(1 + (rand(nSubjects,1) - 0.5)*0.25, weightDifference, 70, 'k', 'filled', ...
            'MarkerFaceAlpha', 0.5);

        yline(0,'k-');
        xlim([0.4 1.6])
        xticks(1)
        xticklabels({'Neural DVs'})
        ylabel({'Difference in absolute values', 'of winning and losing DV weights'})
        title(sprintf('Stimulus offset, p = %.2g', pValue_difference), 'FontWeight', 'normal')
        box off
        set(gca,'FontSize',14)

    end

end
