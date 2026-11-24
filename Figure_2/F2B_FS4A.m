%% =======================================================================
% Figure 2B (PMd) and Figure S4A (M1)
%
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
decoderFolder = fullfile(projectFolder, 'data', 'decoder_output', 'ROI_bilateral_all_areas');

% Path to functions
addpath(fullfile(projectFolder, 'functions'));


%% =======================================================================
% SETTINGS
%% =======================================================================

subjects = { ...
    'S01','S04','S06','S07','S12','S15','S17','S18','S11','S05', ...
    'S19','S20','S22','S10','S02','S08','S24','S21','S23','S25'};

nSubjects     = numel(subjects);
nPermutations = 10000;    % permutations of the permutation tests

% ROIs: name and areas whose decoder output is averaged
roiNames = {'PMd', 'M1'};
roiAreas = {{'6a','6d'}, {'4'}};

% Time: 200 samples, 10 ms each, first sample at -0.5 s
nTime    = 200;
timeSec  = ((1:nTime) - 1) * 0.01 - 0.5;

% Ramp window (s)
rampWindow = timeSec >= 0.2 & timeSec <= 1;

% Evidence-strength bins (|delta contrast|)
evidenceEdges = [0 0.06 0.12 inf];
nEvidenceBins = 3;

% Colors of the evidence bins: low (light orange) to high (dark red)
evidenceColors = [0.99 0.80 0.62; 0.82 0.43 0.34; 0.65 0.05 0.05];




%% =======================================================================
% LOOP OVER ROIs
%% =======================================================================

for roi = 1:numel(roiNames)

    % subjects x time (all trials) and subjects x evidence bins x time
    accuracy             = nan(nSubjects,nTime);
    accuracy_evidenceBin = nan(nSubjects,nEvidenceBins,nTime);

    % Correlation between DV slope and evidence strength
    correlation_slope_evidence = nan(nSubjects,1);


    %% ===================================================================
    % MAIN LOOP ACROSS SUBJECTS
    %% ===================================================================

    for subject = 1:nSubjects

        %% ---------------------------------------------------------------
        % LOAD DECODER OUTPUT (all sessions)
        %
        % logOdds : trials x time, decoder log odds (> 0 = right choice),
        %           averaged across the areas of the ROI
        % choice  : -1 = left, +1 = right
        % stimulus: -1 = left > right, +1 = right > left
        %% ---------------------------------------------------------------

        files = dir(fullfile(decoderFolder, ...
            [subjects{subject} '_*_' roiAreas{roi}{1} '_CHOICEdecoder_alltrials_all.mat']));

        logOdds  = [];
        choice   = [];
        stimulus = [];
        evidence = [];

        for file = 1:numel(files)

            session = sscanf(files(file).name, [subjects{subject} '_%d_']);

            sessionLogOdds = 0;

            for area = 1:numel(roiAreas{roi})
                D = load(fullfile(decoderFolder, sprintf('%s_%d_%s_CHOICEdecoder_alltrials_all.mat', ...
                    subjects{subject}, session, roiAreas{roi}{area})));
                sessionLogOdds = sessionLogOdds + double(D.log_odds_all(:,1:nTime)) / numel(roiAreas{roi});
            end

            logOdds  = [logOdds;  sessionLogOdds];
            choice   = [choice;   double(D.Y_all(:))];
            stimulus = [stimulus; double(D.stimulus_all(:))];
            evidence = [evidence; double(D.evidence_strength(:))];

        end


        %% ---------------------------------------------------------------
        % DECODER ACCURACY
        %% ---------------------------------------------------------------

        isCorrectPrediction = sign(logOdds) == choice;

        accuracy(subject,:) = mean(isCorrectPrediction,1);

        evidenceBin = discretize(abs(evidence), evidenceEdges);

        for bin = 1:nEvidenceBins
            accuracy_evidenceBin(subject,bin,:) = mean(isCorrectPrediction(evidenceBin == bin,:),1);
        end


        %% ---------------------------------------------------------------
        % SINGLE-TRIAL DV SLOPE VS EVIDENCE STRENGTH
        %
        %% ---------------------------------------------------------------

        rampTime   = timeSec(rampWindow) - mean(timeSec(rampWindow));
        trialSlope = (logOdds(:,rampWindow) * rampTime') ./ sum(rampTime.^2);

        correctTrials = stimulus == choice;

        correlation_slope_evidence(subject) = corr( ...
            abs(evidence(correctTrials)), ...
            trialSlope(correctTrials) .* stimulus(correctTrials), ...
            'Type', 'Spearman');

    end


    %% ===================================================================
    % STATISTICS
    %
    % Accuracy against chance (0.5) at each time point, FDR-corrected;
    % slope-evidence correlation against zero
    %% ===================================================================

    pValues = nan(1,nTime);

    for time = 1:nTime
        pValues(time) = permutationTest( ...
            accuracy(:,time) - 0.5, ...
            zeros(nSubjects,1), ...
            nPermutations);
    end

    significantTime = double(fdr_bh(pValues,0.05));
    significantTime(significantTime == 0) = nan;

    pValue_correlation = permutationTest( ...
        correlation_slope_evidence, ...
        zeros(nSubjects,1), ...
        nPermutations);


    %% ===================================================================
    % FIGURE (Figure 2B for PMd, Figure S4A for M1)
    %% ===================================================================

    figure('Color','w','Position',[100 100 600 600],'Name',roiNames{roi});

    % ---------------------------------------------------------------
    % Bottom: accuracy over time (all trials)
    % ---------------------------------------------------------------

    axes('Position',[0.10 0.14 0.60 0.28]);
    hold on

    shadedErrorBar( ...
        timeSec, ...
        nanmean(accuracy,1), ...
        nanstd(accuracy,0,1) ./ sqrt(nSubjects), ...
        'lineProps', {'-','color','k'});

    % Significant time points (p < 0.05, FDR-corrected)
    scatter(timeSec, significantTime .* 0.5, 60, 'k', 'filled', 'square');

    % Stimulus onset, offset and sample onsets
    xline(0, '--', 'LineWidth', 3, 'Color', [0.5 0.5 0.5]);
    xline(1, '--', 'LineWidth', 3, 'Color', [0.5 0.5 0.5]);
    for sample = 1:9
        xline(0.1*sample, '--', 'LineWidth', 2, 'Color', [0.5 0.5 0.5]);
    end

    xlim([-0.5 1.3])
    ylim([0.48 0.8])
    xlabel('Time from stimulus onset (s)')
    ylabel({'Decoder accuracy'; '(fraction correct)'})
    set(gca,'TickDir','out','Box','off','FontSize',15)

  
    axes('Position',[0.24 0.56 0.36 0.36]);
    hold on

    for bin = 1:nEvidenceBins

        binAccuracy = squeeze(nanmean(accuracy_evidenceBin(:,bin,:),1))';

        smoothedAccuracy = smoothdata(binAccuracy, 'gaussian', 10);
        plot(timeSec(rampWindow), smoothedAccuracy(rampWindow), ...
            '--', 'Color', evidenceColors(bin,:), 'LineWidth', 2);

        fitParameters = polyfit(timeSec(rampWindow), binAccuracy(rampWindow), 1);
        plot(timeSec(rampWindow), polyval(fitParameters, timeSec(rampWindow)), ...
            'Color', evidenceColors(bin,:), 'LineWidth', 3);

    end

    xlim([0.2 1])
    title('Evidence strength: low (light) to high (dark)', 'FontWeight', 'normal', 'FontSize', 11)
    set(gca,'TickDir','out','Box','off','FontSize',15)

    % ---------------------------------------------------------------
    % correlation between DV slope and evidence strength
    % ---------------------------------------------------------------

    axes('Position',[0.80 0.56 0.13 0.30]);
    hold on

    bar(1, mean(correlation_slope_evidence), 'FaceColor', 'k', 'FaceAlpha', 0.5);

    errorbar(1, mean(correlation_slope_evidence), ...
        std(correlation_slope_evidence) / sqrt(nSubjects), ...
        'k', 'LineWidth', 1.5);

    scatter(1 + 0.12*(rand(nSubjects,1) - 0.5), correlation_slope_evidence, ...
        60, 'k', 'filled', 'MarkerFaceAlpha', 0.5);

    xticks([])
    xlim([0.4 1.6])
    ylabel('Correlation coefficient')
    title({'Evidence ~ DV slope'; sprintf('p = %.2g', pValue_correlation)}, ...
        'FontWeight', 'normal', 'FontSize', 12)
    set(gca,'TickDir','out','Box','off','FontSize',15)

end
