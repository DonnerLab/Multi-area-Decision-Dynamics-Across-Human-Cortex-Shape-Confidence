%% =======================================================================
% Figure 4A-D and Figure 5A-B
%
% Neural DVs (bilateral choice decoders) of the 22 selected ROIs
%
%
% PMd/M1 DV for the projection = mean DV of areas 4, 6a and 6d
% In the correlation matrix: PMd = mean DV of 6a and 6d, M1 = area 4
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

decoderFolder   = fullfile(projectFolder, 'data', 'decoder_output', 'ROI_bilateral');
functionsFolder = fullfile(projectFolder, 'functions');

addpath(functionsFolder)
addpath(fullfile(functionsFolder, 'gifti-1.6'))


%% =======================================================================
% SETTINGS
%% =======================================================================

subjects = { ...
    'S01','S04','S06','S07','S12','S15','S17','S18','S11','S05', ...
    'S19','S20','S22','S10','S02','S08','S24','S21','S23','S25'};

nSubjects     = numel(subjects);
nPermutations = 10000;    % permutations of the permutation tests

% The 22 ROIs (name in the figure, area name of the decoder files)
roiNames = {'6a(PMd)','6d(PMd)','M1','5mv','6v','OP1','RI','PFcm','LBelt','PSL','7Pm', ...
            '7AL','7Am','7PL','AIP','PFt','PF','PFm','PCV','7m','p47r','8C'};
roiAreas = [{'6a','6d','4'}, roiNames(4:end)];
isMotor  = ismember(roiAreas, {'6a','6d','4'});
nROIs    = numel(roiNames);

% Areas of the PMd/M1 DV that is projected out (Figure 5)
referenceAreas = {'4','6a','6d'};

% ROIs of the correlation matrix: association ROIs, PMd and M1, and their
% order in the figure
matrixNames = [roiNames(~isMotor), {'PMd','M1'}];
matrixAreas = [num2cell(roiAreas(~isMotor)), {{'6a','6d'}}, {{'4'}}];
nMatrix     = numel(matrixNames);

matrixOrder = {'p47r','8C','6v','PMd','M1','PFcm','PF','PFm','PFt','AIP','OP1','RI','LBelt', ...
               'PSL','PCV','5mv','7Pm','7PL','7AL','7Am','7m'};
[~, matrixOrder] = ismember(matrixOrder, matrixNames);

% Time: -500 to 1500 ms (same number of samples as the decoder output)
firstFile = dir(fullfile(decoderFolder, [subjects{1} '_*_6a_CHOICEdecoder_alltrials_all_b.mat']));
firstFile = load(fullfile(decoderFolder, firstFile(1).name), 'log_odds_all');
nTime     = size(firstFile.log_odds_all,2);
time      = linspace(-500, 1500, nTime);
dt        = mean(diff(time));

% End of stimulus viewing (accuracy, projection, correlation matrix)
endWindow = find(time >= 900 & time <= 1000);

% Ramp onset
slopeWindow    = round(100 / dt);                  % 100 ms local slope
[~, normalization] = min(abs(time - 1350));       % DV normalized by its value at 1350 ms
peakWindow     = time >= 1000 & time <= 1400;      % the DV must rise until its peak in this window
nBootstrap     = 1000;                             % bootstrap samples of the baseline slope SD
nSD            = 2;                                % threshold: baseline + 2 SD
minPositive    = 0.75;                             % fraction of positive slopes until the peak

% Least-squares slope of the last slopeWindow samples: slope = values * slopeKernel
centeredTime = ((0:slopeWindow-1)' * dt) - mean((0:slopeWindow-1)' * dt);
slopeKernel  = centeredTime ./ sum(centeredTime.^2);

% Minimum number of trials of a stimulus category x choice (correlations)
minTrials = 10;

% All areas needed for the correlation matrix
neededAreas = unique([referenceAreas, matrixAreas{:}]);

% Reproducible bootstrap and permutations
rng(42)


%% =======================================================================
% ATLAS AND INFLATED LEFT-HEMISPHERE SURFACE
%% =======================================================================

load(fullfile(functionsFolder, 'Glasser_surface.mat'))
atlasTable  = readtable(fullfile(functionsFolder, 'Glasser_labels.csv'));
atlasLabels = atlasTable.Label;

surfaceMesh = gifti(fullfile(functionsFolder, 'S1200.L.inflated_MSMAll.32k_fs_LR.surf.gii'));
faces       = double(surfaceMesh.faces);
vertices    = double(surfaceMesh.vertices);

% Surface vertices of each ROI (left hemisphere = second half of the atlas)
roiVertices = cell(nROIs,1);

for roi = 1:nROIs
    labelIndex = find(strcmp(atlasLabels, ['L_' roiAreas{roi}]));
    roiVertices{roi} = find(atlas.indexmax(32493:end) == labelIndex);
end


%% =======================================================================
% LOOP OVER ANALYSES: 1 = neural DVs (Figure 4), 2 = residual DVs (Figure 5)
%% =======================================================================

for analysis = 1:2

    isResidual = analysis == 2;

    %% ===================================================================
    % DECODING ACCURACY, NORMALIZED DV AND RAMP ONSET OF EACH ROI
    %% ===================================================================

    rampOnset    = nan(nSubjects,nROIs);
    accuracy     = nan(nSubjects,nROIs);
    DV_leftTrials  = nan(nSubjects,nROIs,nTime);   % -DV of left-choice trials
    DV_rightTrials = nan(nSubjects,nROIs,nTime);   %  DV of right-choice trials

    for roi = 1:nROIs

        % Residual DVs: PMd/M1 ROIs not analyzed
        if isResidual && isMotor(roi)
            continue
        end

        for subject = 1:nSubjects

            files = dir(fullfile(decoderFolder, ...
                [subjects{subject} '_*_' roiAreas{roi} '_CHOICEdecoder_alltrials_all_b.mat']));

            signedTrials     = [];   % all trials, sign flipped for left choices
            sessionAccuracy  = [];
            sessionLeft      = [];
            sessionRight     = [];

            for file = 1:numel(files)

                % S02 session 6: decoder output of 7PL, p47r and 8C not
                % reconstructed (MEG data missing for part of the
                % recordings), left out
                if strcmp(subjects{subject}, 'S02') && sscanf(files(file).name, [subjects{subject} '_%d_']) == 6 ...
                        && ismember(roiAreas{roi}, {'7PL','p47r','8C'})
                    continue
                end

                %% -------------------------------------------------------
                % LOAD THIS SESSION
                %
                % log_odds_all: trials x time, DV (> 0 = right choice)
                % Y_all       : -1 = left choice, +1 = right choice
                % accuracy_all: folds x time, cross-validated accuracy
                %% -------------------------------------------------------

                D = load(fullfile(decoderFolder, files(file).name));

                DV     = double(D.log_odds_all);
                choice = double(D.Y_all(:));

                if isResidual

                    % PMd/M1 DV: mean of the reference areas (same session)
                    referenceDV = 0;
                    for area = 1:numel(referenceAreas)
                        R = load(fullfile(decoderFolder, strrep(files(file).name, ...
                            ['_' roiAreas{roi} '_CHOICEdecoder'], ['_' referenceAreas{area} '_CHOICEdecoder'])), 'log_odds_all');
                        referenceDV = referenceDV + double(R.log_odds_all) / numel(referenceAreas);
                    end

                    % Project the PMd/M1 DV out of the ROI DV: weight from
                    % the end of stimulus viewing, applied to all time
                    % points

                    % Normalization: PMd/M1 DV scaled so that its vector
                    % across trials at the end of stimulus viewing has
                    % unit length (reference vector)
                    referenceDV  = referenceDV / norm(mean(referenceDV(:,endWindow),2));
                    referenceEnd = mean(referenceDV(:,endWindow),2);

                    % Projection: ROI DV at the end of stimulus viewing on
                    % the reference vector
                    roiEnd           = mean(DV(:,endWindow),2);
                    projectionWeight = referenceEnd' * roiEnd;

                    % Residual: PMd/M1 DV scaled with the weight removed
                    DV = DV - projectionWeight * referenceDV;

                    % Accuracy of the sign of the residual DV
                    residualEnd = mean(DV(:,endWindow),2);
                    sessionAccuracy(end+1) = mean(sign(residualEnd(residualEnd ~= 0)) == choice(residualEnd ~= 0));

                else

                    % Stored decoder accuracy at the end of stimulus viewing
                    sessionAccuracy(end+1) = nanmean(reshape(D.accuracy_all(:,endWindow),[],1));

                end

                % Mean DV of left (sign flipped) and right choices
                sessionLeft(end+1,:)  = -nanmean(DV(choice == -1,:),1);
                sessionRight(end+1,:) =  nanmean(DV(choice ==  1,:),1);

                % Trials signed towards the choice
                DV(choice == -1,:) = -DV(choice == -1,:);
                signedTrials = [signedTrials; DV];

            end

            if isempty(signedTrials)
                continue
            end

            accuracy(subject,roi)         = mean(sessionAccuracy);
            DV_leftTrials(subject,roi,:)  = nanmean(sessionLeft,1);
            DV_rightTrials(subject,roi,:) = nanmean(sessionRight,1);


            %% -----------------------------------------------------------
            % RAMP ONSET
            %
            % Normalized DV = trial-averaged signed DV (smoothed) divided
            % by its value at 1350 ms. Local slope = least-squares slope
            % in the last 100 ms. Threshold = mean pre-stimulus slope +
            % 2 SD, with SD = bootstrap SD of the pre-stimulus slope
            % (trials resampled). Onset = start of the first 100 ms
            % window (after stimulus onset) whose slope exceeds the
            % threshold, if >= 75% of the slopes from there to the DV
            % peak (1 - 1.4 s) are positive.
            %% -----------------------------------------------------------

            nTrials = size(signedTrials,1);
            meanDV  = smoothdata(nanmean(signedTrials,1), 'gaussian', 20);
            normalizationValue = meanDV(normalization);

            if nTrials < 10 || ~isfinite(normalizationValue) || abs(normalizationValue) <= 1e-6
                continue
            end

            % Bootstrap: trial counts of each bootstrap sample
            % (nBootstrap x trials), mean DV of each bootstrap sample
            bootstrapCounts = histc(randi(nTrials, nTrials, nBootstrap), 1:nTrials)';
            validTrials     = isfinite(signedTrials);
            trialsNoNaN     = signedTrials;
            trialsNoNaN(~validTrials) = 0;

            bootstrapDV = (bootstrapCounts * trialsNoNaN) ./ (bootstrapCounts * double(validTrials));
            bootstrapDV = smoothdata(bootstrapDV, 2, 'gaussian', 20) ./ normalizationValue;

            % Local slopes of the normalized DV and of the bootstrap DVs
            slope          = nan(1,nTime);
            bootstrapSlope = nan(nBootstrap,nTime);

            for t = slopeWindow:nTime
                slopeSamples = (t-slopeWindow+1):t;
                slope(t)            = (meanDV(slopeSamples) ./ normalizationValue) * slopeKernel;
                bootstrapSlope(:,t) = bootstrapDV(:,slopeSamples) * slopeKernel;
            end

            % Threshold
            bootstrapSD = std(bootstrapSlope,0,1);
            bootstrapSD = sqrt(nanmean(bootstrapSD(time <= 0).^2));
            threshold   = nanmean(slope(time <= 0)) + nSD * bootstrapSD;

            % End of the ramp: peak of the DV between 1 and 1.4 s
            peakTimes  = find(peakWindow);
            [~, peak]  = max(meanDV(peakTimes));
            rampEnd    = peakTimes(peak);

            % Start time of the slope window ending at each time point
            windowStart = nan(1,nTime);
            windowStart(slopeWindow:end) = time(1:end-slopeWindow+1);

            crossings = find(slope > threshold & (1:nTime) <= rampEnd & windowStart >= 0);

            for crossing = crossings
                if mean(slope(crossing:rampEnd) > 0) >= minPositive
                    rampOnset(subject,roi) = windowStart(crossing);
                    break
                end
            end

        end

    end

    meanOnset    = nanmean(rampOnset,1);
    semOnset     = nanstd(rampOnset,0,1) ./ sqrt(sum(isfinite(rampOnset),1));
    meanAccuracy = nanmean(accuracy,1);

    % Color of each ROI by rank of ramp onset (early = light blue)
    nWithOnset = sum(isfinite(meanOnset));
    [~, onsetOrder] = sort(meanOnset);
    onsetOrder = onsetOrder(1:nWithOnset);

    blues    = brewermap(256, 'Blues');
    roiColor = repmat([0.6 0.6 0.6], nROIs, 1);
    if nWithOnset > 0
        colorLevels = round(linspace(0.35*256, 256, max(nWithOnset,2)));
        roiColor(onsetOrder,:) = blues(colorLevels(1:nWithOnset),:);
    end


    %% ===================================================================
    % GROUP NORMALIZED DV (Figure 4B)
    %
    % Left and right DVs smoothed and averaged, normalized by the value
    % at 1350 ms, smoothed again for display
    %% ===================================================================

    normalizedDV = nan(nROIs,nTime);

    for roi = 1:nROIs

        leftDV  = smoothdata(squeeze(DV_leftTrials(:,roi,:)),  2, 'gaussian', 20);
        rightDV = smoothdata(squeeze(DV_rightTrials(:,roi,:)), 2, 'gaussian', 20);

        groupDV = (nanmean(leftDV,1) + nanmean(rightDV,1)) ./ 2;

        if isfinite(groupDV(normalization)) && abs(groupDV(normalization)) > 1e-6
            normalizedDV(roi,:) = smoothdata(groupDV ./ groupDV(normalization), 2, 'gaussian', 10);
        end

    end


    %% ===================================================================
    % CORRELATION MATRIX (Figure 4D, 5B)
    %
    % Per session: DV of each ROI averaged over 0.9 - 1 s; correlation
    % across trials within each stimulus category x choice (>= 10
    % trials), averaged across categories, then sessions, then subjects
    %% ===================================================================

    % subjects x ROI x ROI
    correlationMatrix = nan(nSubjects,nMatrix,nMatrix);

    for subject = 1:nSubjects

        files = dir(fullfile(decoderFolder, [subjects{subject} '_*_6a_CHOICEdecoder_alltrials_all_b.mat']));

        sessionMatrix = nan(numel(files),nMatrix,nMatrix);

        for file = 1:numel(files)

            D        = load(fullfile(decoderFolder, files(file).name), 'Y_all', 'stimulus_all');
            choice   = double(D.Y_all(:));
            stimulus = sign(double(D.stimulus_all(:)));

            % DV of every area of this session at the end of stimulus
            % viewing (trials x areas)
            endDV = nan(numel(choice),numel(neededAreas));

            for area = 1:numel(neededAreas)

                % S02 session 6: decoder output of 7PL, p47r and 8C not
                % reconstructed (MEG data missing for part of the
                % recordings), left out (NaN)
                if strcmp(subjects{subject}, 'S02') && sscanf(files(file).name, [subjects{subject} '_%d_']) == 6 ...
                        && ismember(neededAreas{area}, {'7PL','p47r','8C'})
                    continue
                end

                A = load(fullfile(decoderFolder, strrep(files(file).name, ...
                    '_6a_CHOICEdecoder', ['_' neededAreas{area} '_CHOICEdecoder'])), 'log_odds_all');
                endDV(:,area) = mean(double(A.log_odds_all(:,endWindow)),2);
            end

            % PMd/M1 DV (residual analysis)
            if isResidual
                referenceDV = mean(endDV(:,ismember(neededAreas,referenceAreas)),2);

                % Normalization: PMd/M1 DV across the trials of the
                % session scaled to unit length (reference vector)
                reference = referenceDV / norm(referenceDV);
            end

            % DV of each ROI of the matrix (mean of its areas)
            roiDV = nan(numel(choice),nMatrix);

            for roi = 1:nMatrix

                if isResidual && any(strcmp(matrixNames{roi}, {'PMd','M1'}))
                    continue
                end

                areaDV = mean(endDV(:,ismember(neededAreas,matrixAreas{roi})),2);

                % Residual: PMd/M1 DV projected out
                if isResidual

                    % Projection: ROI DV on the reference vector
                    projectionWeight = reference' * areaDV;

                    % Residual: reference vector scaled with the weight removed
                    areaDV = areaDV - projectionWeight * reference;

                end

                roiDV(:,roi) = areaDV;

            end

            % Correlations within each stimulus category x choice
            categoryMatrix = nan(nMatrix,nMatrix,4);
            category = 0;

            for choiceSide = [-1 1]
                for stimulusSide = [-1 1]

                    category = category + 1;
                    trials = choice == choiceSide & stimulus == stimulusSide;

                    if sum(trials) >= minTrials
                        categoryMatrix(:,:,category) = corr(roiDV(trials,:), 'Rows', 'pairwise');
                    end

                end
            end

            sessionMatrix(file,:,:) = nanmean(categoryMatrix,3);

        end

        correlationMatrix(subject,:,:) = nanmean(sessionMatrix,1);

    end

    % Group average and two-sided permutation tests against zero,
    % FDR-corrected across ROI pairs
    groupMatrix = squeeze(nanmean(correlationMatrix,1));

    upperTriangle = find(triu(true(nMatrix),1) & isfinite(groupMatrix));
    pValues = nan(size(upperTriangle));

    for pair = 1:numel(upperTriangle)
        [roi1, roi2] = ind2sub([nMatrix nMatrix], upperTriangle(pair));
        pairCorrelation = correlationMatrix(:,roi1,roi2);
        pairCorrelation = pairCorrelation(isfinite(pairCorrelation));
        pValues(pair) = permutationTest(pairCorrelation, zeros(size(pairCorrelation)), nPermutations);
    end

    significantPairs = false(nMatrix);
    significantPairs(upperTriangle) = fdr_bh(pValues, 0.05);
    significantPairs = significantPairs | significantPairs';


    %% ===================================================================
    % FIGURES
    %% ===================================================================

    % Panels: A = map (lateral, medial), B = normalized DV, C = ramp
    % onset, D = correlation matrix (Figure 5: A = map, B = matrix)
    if isResidual
        figure('Color','w','Name','Figure 5A-B','Position',[40 40 1000 450]);
        nRows = 2;  mapPanels = [1 3];  matrixPanel = [2 4];
    else
        figure('Color','w','Name','Figure 4A-D','Position',[40 40 1000 900]);
        nRows = 4;  mapPanels = [1 3];  matrixPanel = [6 8];
    end
    nColumns = 2;

    % ---------------------------------------------------------------
    % A: decoding accuracy (lateral and medial view), ROI labels in the
    %    color of their ramp onset
    % ---------------------------------------------------------------

    accuracyColors = brewermap(256, 'Reds');
    accuracyColors = accuracyColors(round(0.15*256):end,:);
    accuracyLimits = [0.5 0.55];

    vertexColors = repmat([0.7 0.7 0.7], size(vertices,1), 1);

    for roi = 1:nROIs
        if isfinite(meanAccuracy(roi)) && ~isempty(roiVertices{roi})
            colorIndex = round(1 + (meanAccuracy(roi) - accuracyLimits(1)) / diff(accuracyLimits) * (size(accuracyColors,1)-1));
            colorIndex = min(max(colorIndex,1), size(accuracyColors,1));
            vertexColors(roiVertices{roi},:) = repmat(accuracyColors(colorIndex,:), numel(roiVertices{roi}), 1);
        end
    end

    viewAngles = [270 0; 90 0];   % lateral, medial

    for view_n = 1:2

        subplot(nRows,nColumns,mapPanels(view_n))

        patch('Faces', faces, 'Vertices', vertices, 'FaceVertexCData', vertexColors, ...
            'FaceColor', 'interp', 'EdgeColor', 'none');

        axis equal
        axis off
        view(viewAngles(view_n,:))
        lighting flat
        material([0.85 0.10 0 1 1])
        camlight('headlight')

        % Labels of the ROIs seen from this side, at the ROI center
        for roi = 1:nROIs
            if isfinite(meanAccuracy(roi)) && ~isempty(roiVertices{roi})
                roiCenter = mean(vertices(roiVertices{roi},:),1);
                isLateral = roiCenter(1) < mean(vertices(:,1));
                if isLateral == (view_n == 1)
                    text(roiCenter(1) + 40*(2*(view_n == 2) - 1), roiCenter(2), roiCenter(3), roiNames{roi}, ...
                        'Color', roiColor(roi,:), 'FontWeight', 'bold', 'FontSize', 9, ...
                        'HorizontalAlignment', 'center');
                end
            end
        end

    end

    colormap(gca, accuracyColors);
    caxis(accuracyLimits);
    colorbar;

    if ~isResidual

        % -----------------------------------------------------------
        % B: normalized DV of each ROI (colored by ramp onset)
        % -----------------------------------------------------------

        subplot(nRows,nColumns,[2 4])
        hold on

        for roi = 1:nROIs
            plot(time, normalizedDV(roi,:), 'Color', roiColor(roi,:), 'LineWidth', 1.2);
        end

        xline(0, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.5);
        xline(1000, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.5);
        yline(1, ':', 'Color', [0.4 0.4 0.4]);
        xlim([-500 1400])
        xticks([-500 0 500 1000])
        xlabel('Time from stimulus onset (ms)')
        ylabel('Normalized decision variable')
        set(gca,'TickDir','out','Box','off','FontSize',10)

        % -----------------------------------------------------------
        % C: ramp onset (mean +/- SEM), ROIs sorted by onset
        % -----------------------------------------------------------

        subplot(nRows,nColumns,[5 7])
        hold on

        for rank = 1:nWithOnset
            roi = onsetOrder(rank);
            errorbar(meanOnset(roi), rank, semOnset(roi), 'horizontal', 'o', ...
                'Color', roiColor(roi,:), 'MarkerFaceColor', roiColor(roi,:), ...
                'MarkerEdgeColor', 'k', 'LineWidth', 1.2);
        end

        set(gca, 'YDir', 'reverse', 'YTick', 1:nWithOnset, 'YTickLabel', roiNames(onsetOrder), ...
            'TickLabelInterpreter', 'none', 'TickDir', 'out', 'Box', 'off', 'FontSize', 9);
        xlabel('Ramp onset (ms)')

    end

    % ---------------------------------------------------------------
    % D (Figure 4) / B (Figure 5): correlation matrix; non-significant
    % pairs transparent, PMd and M1 rows outlined
    % ---------------------------------------------------------------

    subplot(nRows,nColumns,matrixPanel)
    hold on

    orderedMatrix      = groupMatrix(matrixOrder,matrixOrder);
    orderedSignificant = significantPairs(matrixOrder,matrixOrder) | eye(nMatrix);

    transparency = 0.35 * isfinite(orderedMatrix);
    transparency(orderedSignificant & isfinite(orderedMatrix)) = 1;

    orderedMatrix(~isfinite(orderedMatrix)) = 0;

    imagesc(orderedMatrix, 'AlphaData', transparency);
    axis image
    set(gca, 'YDir', 'reverse')
    colormap(gca, flipud(brewermap(256, 'RdBu')));
    caxis([-1 1]);
    colorbar;

    for row = find(ismember(matrixNames(matrixOrder), {'PMd','M1'}))
        rectangle('Position', [0.5 row-0.5 nMatrix 1], 'EdgeColor', 'k', 'LineWidth', 1.5);
    end

    set(gca, 'XTick', 1:nMatrix, 'XTickLabel', matrixNames(matrixOrder), ...
        'YTick', 1:nMatrix, 'YTickLabel', matrixNames(matrixOrder), ...
        'TickLabelInterpreter', 'none', 'TickDir', 'out', 'FontSize', 8);
    xtickangle(90)


   
end