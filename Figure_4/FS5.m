%% =======================================================================
% Figure S5
%
% Neural DVs of the selected ROIs (bilateral choice decoders)
%
% S5A: correlation between the slope of the neural DV during stimulus
%      viewing (0.2 - 1 s) and evidence strength (criterion 2 of
%      ROI_selection_TableS1.m); bars = group mean +/- SEM
% S5B: across-trial variance of the neural DV over time (group mean +/-
%      SEM); red = straight line fitted to the group mean in the second
%      half of stimulus viewing (0.5 - 1 s)
%
% PMd = mean DV of areas 6a and 6d, M1 = area 4
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

% Bilateral decoder output of the 22 ROIs (same files as Figures 4, 5 and Table S1)
decoderFolder = fullfile(projectFolder, 'data', 'decoder_output', 'ROI_bilateral');

% Path to functions
addpath(fullfile(projectFolder, 'functions'));


%% =======================================================================
% SETTINGS
%% =======================================================================

subjects = { ...
    'S01','S04','S06','S07','S12','S15','S17','S18','S11','S05', ...
    'S19','S20','S22','S10','S02','S08','S24','S21','S23','S25'};

nSubjects = numel(subjects);

% ROIs (name, areas whose DVs are averaged, cluster), in the order of the
% correlation matrix of Figure 4D
roiNames = {'p47r','8C','6v','PMd','M1','PFcm','PF','PFm','PFt','AIP','OP1','RI','LBelt', ...
            'PSL','PCV','5mv','7Pm','7PL','7m','7AL','7Am'};
roiAreas = num2cell(roiNames);
roiAreas{strcmp(roiNames,'PMd')} = {'6a','6d'};
roiAreas{strcmp(roiNames,'M1')}  = {'4'};
roiAreas = cellfun(@cellstr, roiAreas, 'UniformOutput', false);
nROIs    = numel(roiNames);

% Cluster colors: prefrontal (orange), PMd/M1 (black), lateral parietal
% (blue), medial parietal (green)
roiCluster    = [1 1 1 2 2 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4];
clusterColors = [0.93 0.49 0.13; 0 0 0; 0.25 0.58 0.85; 0.27 0.66 0.30];
roiColors     = clusterColors(roiCluster,:);

% Time: 200 samples, 10 ms each, first sample at -500 ms
timeMs = ((1:200) - 1) * 10 - 500;

slopeWindow    = 71:151;                             % 200 - 1000 ms
varianceWindow = find(timeMs >= 500 & timeMs <= 1000);

% ROIs x subjects (x time)
evidenceScaling  = nan(nROIs,nSubjects);
varianceOverTime = nan(nROIs,nSubjects,200);


%% =======================================================================
% LOOP OVER ROIs AND SUBJECTS
%% =======================================================================

for roi = 1:nROIs
    for subject = 1:nSubjects

        files = dir(fullfile(decoderFolder, ...
            [subjects{subject} '_*_' roiAreas{roi}{1} '_CHOICEdecoder_alltrials_all_b.mat']));

        sessionVariance = nan(numel(files),200);

        % Pooled z-scored signed evidence and DV slopes of correct trials,
        % left (1) and right (2) choices
        pooledEvidence = {[], []};
        pooledSlope    = {[], []};

        for file = 1:numel(files)

            % DV of this session, averaged across the areas of the ROI
            session = sscanf(files(file).name, [subjects{subject} '_%d_']);

            % S02 session 6: decoder output of 7PL, p47r and 8C not
            % reconstructed (MEG data missing for part of the recordings),
            % left out
            if strcmp(subjects{subject}, 'S02') && session == 6 && ismember(roiNames{roi}, {'7PL','p47r','8C'})
                continue
            end

            DV = 0;
            for area = 1:numel(roiAreas{roi})
                D = load(fullfile(decoderFolder, sprintf('%s_%d_%s_CHOICEdecoder_alltrials_all_b.mat', ...
                    subjects{subject}, session, roiAreas{roi}{area})));
                DV = DV + double(D.log_odds_all(:,1:200)) / numel(roiAreas{roi});
            end

            choice   = sign(double(D.Y_all(:)));
            stimulus = double(D.stimulus_all(:));
            evidence = double(D.evidence_strength(:));

            % Variance across trials at each time point (Figure S5B)
            sessionVariance(file,:) = var(DV,0,1);

            % Slope of the DV of each trial (Figure S5A)
            trialSlope = nan(size(DV,1),1);
            for trial = 1:size(DV,1)
                fitParameters = polyfit(timeMs(slopeWindow), DV(trial,slopeWindow), 1);
                trialSlope(trial) = fitParameters(1);
            end

            signedEvidence = evidence .* stimulus;

            for choiceSide = [-1 1]

                trials = choice == choiceSide & stimulus == choiceSide;

                if sum(trials) >= 6 && std(signedEvidence(trials)) > 0 && std(trialSlope(trials)) > 0
                    side = (choiceSide+3)/2;
                    pooledEvidence{side} = [pooledEvidence{side}; zscore(signedEvidence(trials))];
                    pooledSlope{side}    = [pooledSlope{side};    zscore(trialSlope(trials))];
                end

            end

        end

        % Subject values: mean of left and right correlation; mean
        % variance across sessions
        sideCorrelation = nan(1,2);
        for side = 1:2
            if numel(pooledEvidence{side}) >= 10
                sideCorrelation(side) = corr(pooledEvidence{side}, pooledSlope{side});
            end
        end

        evidenceScaling(roi,subject)    = nanmean(sideCorrelation);
        varianceOverTime(roi,subject,:) = nanmean(sessionVariance,1);

    end
end


%% =======================================================================
% FIGURE 1 (Figure S5A) - DV SLOPE VS EVIDENCE
%% =======================================================================

figure(1)
hold on

for roi = 1:nROIs

    bar(roi, nanmean(evidenceScaling(roi,:)), 0.7, 'FaceColor', roiColors(roi,:), 'EdgeColor', 'k');

    errorbar(roi, nanmean(evidenceScaling(roi,:)), ...
        nanstd(evidenceScaling(roi,:)) / sqrt(sum(isfinite(evidenceScaling(roi,:)))), ...
        'Color', [0.25 0.25 0.25], 'LineStyle', 'none', 'LineWidth', 1.1);

end

yline(0,'k-');
set(gca, 'XTick', 1:nROIs, 'XTickLabel', roiNames, 'XTickLabelRotation', 45, ...
    'FontSize', 10, 'Box', 'off', 'TickDir', 'out');
xlim([0.4 nROIs+0.6])
ylabel('Correlation: DV slope vs evidence (r)')


%% =======================================================================
% FIGURE 2 (Figure S5B) - VARIANCE OF THE NEURAL DV, ONE PANEL PER ROI
%% =======================================================================

figure(2)

for roi = 1:nROIs

    subplot(4,6,roi)
    hold on

    roiVariance = squeeze(varianceOverTime(roi,:,:));    % subjects x time
    meanVariance = nanmean(roiVariance,1);

    shadedErrorBar(timeMs, meanVariance, ...
        nanstd(roiVariance,0,1) ./ sqrt(sum(isfinite(roiVariance),1)), ...
        'lineProps', {'-','color','k'});

    % Straight line fitted to the group mean, 500 - 1000 ms
    fitParameters = polyfit(timeMs(varianceWindow), meanVariance(varianceWindow), 1);
    plot(timeMs(varianceWindow), polyval(fitParameters, timeMs(varianceWindow)), ...
        '-', 'Color', [0.80 0.10 0.10], 'LineWidth', 1.2);

    xline(0, ':', 'Color', [0.5 0.5 0.5]);
    xline(1000, ':', 'Color', [0.5 0.5 0.5]);

    xlim([-500 1490])
    xticks([-500 0 500 1000])
    title(roiNames{roi}, 'Color', roiColors(roi,:))
    box off
    set(gca, 'FontSize', 7, 'TickDir', 'out')

end

subplot(4,6,19)
xlabel('Time from stimulus onset (ms)')
ylabel('Variance of neural DV')