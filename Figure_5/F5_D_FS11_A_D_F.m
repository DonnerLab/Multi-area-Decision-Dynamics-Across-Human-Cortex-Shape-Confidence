%% =======================================================================
% Figure 5D and Figure S11A-D, F
%
% Residual neural DVs of the 19 association ROIs (bilateral choice
% decoders) after projecting out the PMd/M1 DV (mean DV of areas 4, 6a
% and 6d), in sliding windows (350 ms, steps of 100 ms)
%
% Projection, within each session (regression without intercept):
%   residual DV = DV - (PMd/M1' DV / PMd/M1' PMd/M1) PMd/M1
%
% Winning areas: residual DV with the sign of the choice; losing areas:
% opposite sign. Strength of winning (losing) areas = mean |residual DV|
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

% Bilateral decoder output:
decoderFolder = fullfile(projectFolder, 'data', 'decoder_output', 'ROI_bilateral');
% Path to functions
addpath(fullfile(projectFolder, 'functions'));


%% =======================================================================
% SETTINGS
%% =======================================================================

subjects = { ...
    'S01','S04','S06','S07','S12','S15','S17','S18','S11','S05', ...
    'S19','S20','S22','S10','S02','S08','S24','S21','S23','S25'};

nSubjects     = numel(subjects);
nPermutations = 10000;    % permutations of the cluster-based tests

% Association ROIs (selected by ROI_selection_TableS1.m) and PMd/M1 areas
areas          = {'5mv','6v','OP1','RI','PFcm','LBelt','PSL','7Pm','7AL','7Am', ...
                  '7PL','AIP','PFt','PF','PFm','PCV','7m','p47r','8C'};
referenceAreas = {'4','6a','6d'};
nAreas         = numel(areas);

% Sliding windows of 35 samples (350 ms), steps of 10 samples (100 ms);
% decoder time base: 10 ms samples, first sample at -0.5 s
windowLength = 35;
windowStarts = 1:10:(200 - windowLength + 1);
nWindows     = numel(windowStarts);
windowTime   = ((windowStarts + windowLength/2 - 1) * 10 - 500) / 1000;



% subjects x windows (x areas)
correlation_withPMd_before = nan(nSubjects,nWindows,nAreas);   % S11A
correlation_withPMd_after  = nan(nSubjects,nWindows,nAreas);   % S11B
correlation_winLose_before = nan(nSubjects,nWindows);          % S11C, black
correlation_winLose_after  = nan(nSubjects,nWindows);          % S11C, red
correlation_winLose_abs    = nan(nSubjects,nWindows);          % S11D

% subjects x windows x [evidence, winning areas, losing areas]
betas_confidence = nan(nSubjects,nWindows,3);                  % 5D
betas_RT         = nan(nSubjects,nWindows,3);                  % S11F


%% =======================================================================
% MAIN LOOP ACROSS SUBJECTS
%% =======================================================================

for subject = 1:nSubjects

    %% -------------------------------------------------------------------
    % LOAD ALL SESSIONS
    %
    % areaDV{area}: trials x time; referenceDV: mean of the PMd/M1 areas
    % Y_all: -1 = left, +1 = right choice; confidence_all: 1 = low,
    % 2 = high; rts_all: RT; evidence_strength
    %% -------------------------------------------------------------------

    files = dir(fullfile(decoderFolder, [subjects{subject} '_*_' areas{1} '_CHOICEdecoder_alltrials_all_b.mat']));

    areaDV      = cell(1,nAreas);
    referenceDV = [];
    choice      = [];
    confidence  = [];
    RT          = [];
    evidence    = [];
    session     = [];

    for file = 1:numel(files)

        sessionNumber = sscanf(files(file).name, [subjects{subject} '_%d_']);

        % S02 session 6 left out: decoder output of 7PL, p47r and 8C not
        % reconstructed (MEG data missing for part of the recordings)
        if strcmp(subjects{subject}, 'S02') && sessionNumber == 6
            continue
        end

        for area = 1:nAreas
            D = load(fullfile(decoderFolder, sprintf('%s_%d_%s_CHOICEdecoder_alltrials_all_b.mat', ...
                subjects{subject}, sessionNumber, areas{area})));
            areaDV{area} = [areaDV{area}; double(D.log_odds_all(:,1:200))];
        end

        sessionReference = 0;
        for area = 1:numel(referenceAreas)
            R = load(fullfile(decoderFolder, sprintf('%s_%d_%s_CHOICEdecoder_alltrials_all_b.mat', ...
                subjects{subject}, sessionNumber, referenceAreas{area})), 'log_odds_all');
            sessionReference = sessionReference + double(R.log_odds_all(:,1:200)) / numel(referenceAreas);
        end

        referenceDV = [referenceDV; sessionReference];
        choice      = [choice;      sign(double(D.Y_all(:)))];
        confidence  = [confidence;  double(D.confidence_all(:)) - 1];
        RT          = [RT;          double(D.rts_all(:))];
        evidence    = [evidence;    double(D.evidence_strength(:))];
        session     = [session;     file * ones(numel(D.Y_all),1)];

    end

    nTrials = numel(choice);


    %% -------------------------------------------------------------------
    % LOOP OVER WINDOWS
    %% -------------------------------------------------------------------

    for window = 1:nWindows

        windowSamples = windowStarts(window):(windowStarts(window) + windowLength - 1);

        % Window-averaged DVs: trials x areas, and PMd/M1 DV
        DV = nan(nTrials,nAreas);
        for area = 1:nAreas
            DV(:,area) = mean(areaDV{area}(:,windowSamples),2);
        end

        pmdDV = mean(referenceDV(:,windowSamples),2);

        % Residual DVs: PMd/M1 DV projected out within each session
        residualDV = DV;

       for thisSession = unique(session)'

    trials = session == thisSession;

    % Normalization: PMd/M1 DV across the trials of the session,
    % scaled to unit length (reference vector)
    reference = pmdDV(trials) / norm(pmdDV(trials));

    % Projection: one weight per area (area DV on the reference)
    projectionWeight = reference' * DV(trials,:);

    % Residual: reference vector scaled with the weight removed
    residualDV(trials,:) = DV(trials,:) - reference * projectionWeight;

end


        %% ---------------------------------------------------------------
        % S11A-B: CORRELATION OF EACH ROI WITH THE PMd/M1 DV
        %% ---------------------------------------------------------------

        for area = 1:nAreas
            correlation_withPMd_before(subject,window,area) = corr(pmdDV, DV(:,area));
            correlation_withPMd_after(subject,window,area)  = corr(pmdDV, residualDV(:,area));
        end


        %% ---------------------------------------------------------------
        % WINNING AND LOSING AREAS
        %
        % before: signed DVs; after: signed and rectified residual DVs
        %% ---------------------------------------------------------------

        isWinning = sign(DV) == choice;
        isLosing  = sign(DV) ~= choice & DV ~= 0;

        winning_before = DV;  winning_before(~isWinning) = NaN;
        losing_before  = DV;  losing_before(~isLosing)   = NaN;

        isWinning = sign(residualDV) == choice;
        isLosing  = sign(residualDV) ~= choice & residualDV ~= 0;

        winning_after = residualDV;       winning_after(~isWinning) = NaN;
        losing_after  = residualDV;       losing_after(~isLosing)   = NaN;
        winning_abs   = abs(residualDV);  winning_abs(~isWinning)   = NaN;
        losing_abs    = abs(residualDV);  losing_abs(~isLosing)     = NaN;

        % Mean over the winning / losing areas of each trial
        winLose = [nanmean(winning_before,2), nanmean(losing_before,2), ...
                   nanmean(winning_after,2),  nanmean(losing_after,2), ...
                   nanmean(winning_abs,2),    nanmean(losing_abs,2)];


        %% ---------------------------------------------------------------
        % Z-SCORING WITHIN SESSION
     %%------------------------------------------------------------

        winLose_z    = nan(size(winLose));
        regressors_z = nan(nTrials,4);   % evidence, winning, losing, RT

        for thisSession = unique(session)'

            trials = session == thisSession;

            for pair = 1:3
                valid = trials & all(isfinite(winLose(:,2*pair-1:2*pair)),2);
                if sum(valid) > 2
                    winLose_z(valid,2*pair-1:2*pair) = zscore(winLose(valid,2*pair-1:2*pair));
                end
            end

            valid = trials & all(isfinite([winLose(:,5:6) evidence RT confidence]),2);
            if sum(valid) > 2
                regressors_z(valid,:) = zscore([evidence(valid) winLose(valid,5:6) RT(valid)]);
            end

        end

        % S11C-D: correlation between winning and losing areas
        valid = all(isfinite(winLose_z(:,1:2)),2);
        correlation_winLose_before(subject,window) = corr(winLose_z(valid,1), winLose_z(valid,2));

        valid = all(isfinite(winLose_z(:,3:4)),2);
        correlation_winLose_after(subject,window) = corr(winLose_z(valid,3), winLose_z(valid,4));

        valid = all(isfinite(winLose_z(:,5:6)),2);
        correlation_winLose_abs(subject,window) = corr(winLose_z(valid,5), winLose_z(valid,6));


        %% ---------------------------------------------------------------
        % REGRESSIONS (Figure 5D, S11F)
        %% ---------------------------------------------------------------

        valid = RT >= 0 & all(isfinite(regressors_z),2) & isfinite(confidence);

        model_confidence = fitglm(regressors_z(valid,1:3), confidence(valid), 'Distribution', 'binomial');
        betas_confidence(subject,window,:) = model_confidence.Coefficients.Estimate(2:4);

        model_RT = fitlm(regressors_z(valid,1:3), regressors_z(valid,4));
        betas_RT(subject,window,:) = model_RT.Coefficients.Estimate(2:4);

    end

end


%% =======================================================================
% FIGURE 1 (Figure S11A-D)
%% =======================================================================

areaColors = brewermap(nAreas, '*RdBu');

figure('Color','w','Position',[40 40 1300 800]);

% ---------------------------------------------------------------
% A-B: correlation of each ROI with the PMd/M1 DV
% ---------------------------------------------------------------

panelData   = {correlation_withPMd_before, correlation_withPMd_after};
panelTitles = {'A: before removal of PMd/M1 DV', 'B: after removal of PMd/M1 DV'};

for panel = 1:2

    subplot(2,2,panel)
    hold on

    for area = 1:nAreas
        shadedErrorBar(windowTime, nanmean(panelData{panel}(:,:,area),1), ...
            nanstd(panelData{panel}(:,:,area),0,1) ./ sqrt(nSubjects), ...
            'lineProps', {'-','color',areaColors(area,:),'LineWidth',2});
    end

    yline(0,'k--');
    xline(0,'k--');
    xline(1,'k--');
    ylim([-0.05 1])
    xlabel('Time from stimulus onset (s)')
    ylabel('Correlation with PMd/M1 DV')
    title(panelTitles{panel}, 'FontWeight', 'normal')
    box off
    set(gca,'FontSize',11,'TickDir','out')

end

% ---------------------------------------------------------------
% C: winning vs losing areas, before (black) and after (red)
% ---------------------------------------------------------------

subplot(2,2,3)
hold on

shadedErrorBar(windowTime, nanmean(correlation_winLose_before,1), ...
    nanstd(correlation_winLose_before,0,1) ./ sqrt(nSubjects), ...
    'lineProps', {'-','color','k','LineWidth',2});

shadedErrorBar(windowTime, nanmean(correlation_winLose_after,1), ...
    nanstd(correlation_winLose_after,0,1) ./ sqrt(nSubjects), ...
    'lineProps', {'-','color','r','LineWidth',2});

yline(0,'k--');
xline(0,'k--');
xline(1,'k--');
ylim([-1 0])
xlabel('Time from stimulus onset (s)')
ylabel('Correlation coeff.')
title('C: winning vs losing areas, before (black) / after (red) removal', 'FontWeight', 'normal')
box off
set(gca,'FontSize',11,'TickDir','out')

% ---------------------------------------------------------------
% D: rectified residual DVs of winning vs losing areas
% ---------------------------------------------------------------

subplot(2,2,4)
hold on

shadedErrorBar(windowTime, nanmean(correlation_winLose_abs,1), ...
    nanstd(correlation_winLose_abs,0,1) ./ sqrt(nSubjects), ...
    'lineProps', {'-','color','b','LineWidth',2});

yline(0,'k--');
xline(0,'k--');
xline(1,'k--');
ylim([-0.2 0.15])
xlabel('Time from stimulus onset (s)')
ylabel('Correlation coeff.')
title('D: rectified residual DVs of winning vs losing areas', 'FontWeight', 'normal')
box off
set(gca,'FontSize',11,'TickDir','out')


%% =======================================================================
% FIGURE 2 (Figure 5D) AND FIGURE 3 (Figure S11F)
%
% Weights of the winning (black) and losing (blue) areas over time;
% marks: p < 0.05, cluster-based permutation test against zero
%% =======================================================================

allBetas    = {betas_confidence, betas_RT};
figureNames = {'Figure 5D - Confidence', 'Figure S11F - Reaction time'};
areaTypeColors = [0 0 0; 0 0.45 0.90];

for behavior = 1:2

    figure('Color','w','Position',[100 100 900 350],'Name',figureNames{behavior});
    hold on

    for areaType = 1:2

        % windows x subjects (subjects with all windows fitted)
        weights = squeeze(allBetas{behavior}(:,:,areaType+1))';
        weights = weights(:, all(isfinite(weights),1));

        shadedErrorBar(windowTime, mean(weights,2), std(weights,0,2) ./ sqrt(size(weights,2)), ...
            'lineProps', {'-','color',areaTypeColors(areaType,:),'LineWidth',2});

        significantWindows = double(cluster_perm_1samp(weights, nPermutations, 0.05));
        significantWindows(significantWindows == 0) = nan;

        meanWeights  = squeeze(nanmean(allBetas{behavior}(:,:,2:3),1));
        markerHeight = max(abs(meanWeights(:))) * (1.6 - 0.1*areaType);

        plot(windowTime, significantWindows * markerHeight, '-', ...
            'Color', areaTypeColors(areaType,:), 'LineWidth', 6);

    end

    yline(0,'k--');
    xline(0,'k--');
    xline(1,'k--');
    xlabel('Time from stimulus onset (s)')
    ylabel('Weight (a.u.)')
    title([figureNames{behavior} ': winning (black) and losing (blue) areas'], 'FontWeight', 'normal')
    box off
    set(gca,'FontSize',14)

end