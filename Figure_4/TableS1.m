%% =======================================================================
% Table S1
%
% Three measures of the bilateral neural DV of
% the selected 22 ROIs (PMd/M1 areas 6a, 6d, 4 and the 19 association
% ROIs). 
% (1) Variance of the neural DV increases over time: across-trial
%     variance at each time point (0.5 - 1 s), separately for left and
%     right choices; slope of a straight line fitted to the variance;
%     mean of left and right, then of the sessions. Test: slope > 0.
%
% (2) Build-up of the neural DV scales with evidence strength: slope of
%     the DV of each trial (0.2 - 1 s); correct trials, separately for
%     left and right choices; signed evidence and slope z-scored within
%     session, pooled across sessions; Pearson correlation; mean of left
%     and right. Test: correlation > 0.
%
% (3) Residual DV predicts choice: DV at the end of stimulus viewing
%     (0.9 - 1 s) after projecting out the PMd/M1 DV (mean of areas 6a,
%     6d and 4) within each session; accuracy of its sign.
%     Test: accuracy > 0.5.
%
% Table S1 (tableS1): mean and 95% confidence interval across
% participants (mean +/- t(0.975, n-1) x SEM) of the three measures,
% for the PMd/M1 reference areas and the 19 association ROIs
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

% Bilateral decoder output of all areas: <subject>_<session>_<area>_CHOICEdecoder_alltrials_all_b.mat
decoderFolder = fullfile(projectFolder, 'data', 'decoder_output', 'ROI_bilateral');
fileEnding    = '_CHOICEdecoder_alltrials_all_b.mat';

% Path to functions
addpath(fullfile(projectFolder, 'functions'));


%% =======================================================================
% SETTINGS
%% =======================================================================

subjects = { ...
    'S01','S04','S06','S07','S12','S15','S17','S18','S11','S05', ...
    'S19','S20','S22','S10','S02','S08','S24','S21','S23','S25'};

nSubjects     = numel(subjects);

% Areas: reference areas and the 19 association ROIs (rows of Table S1)
areas = {'6a','6d','4','p47r','8C','6v','PFcm','PFt','PF','PFm','AIP','OP1','RI', ...
         'LBelt','PSL','PCV','5mv','7Pm','7m','7Am','7PL','7AL'};

nAreas = numel(areas);

% Areas of the PMd/M1 DV
referenceAreas = {'6a','6d','4'};

% Rows of Table S1
tableAreas = areas;

% Time: 200 samples, 10 ms each, first sample at -0.5 s
timeSec = ((1:200) - 1) * 0.01 - 0.5;

varianceWindow = find(timeSec >= 0.5 & timeSec <= 1);   % criterion 1
slopeWindow    = 71:151;                                % criterion 2 (0.2 - 1 s)
endWindow      = 141:150;                               % criterion 3 (0.9 - 0.99 s)

% areas x subjects
varianceSlope    = nan(nAreas,nSubjects);
evidenceScaling  = nan(nAreas,nSubjects);
residualAccuracy = nan(nAreas,nSubjects);


%% =======================================================================
% MAIN LOOP ACROSS SUBJECTS
%% =======================================================================

for subject = 1:nSubjects

    % Session numbers of the files of every area of this subject
    areaSessions = cell(1,nAreas);
    for area = 1:nAreas
        files = dir(fullfile(decoderFolder, [subjects{subject} '_*_' areas{area} fileEnding]));
        areaSessions{area} = cellfun(@(name) sscanf(name, [subjects{subject} '_%d_']), {files.name});
    end

    % Sessions with a file of every area (session 1 is saved as 1 or as
    % 11 in the file names: counted as 1)
    sessions = 1:100;
    for area = 1:nAreas
        areaSession = areaSessions{area};
        areaSession(areaSession == 11) = 1;
        sessions = intersect(sessions, areaSession);
    end

    % PMd/M1 DV at the end of stimulus viewing, one cell per session
    referenceDV = cell(1,numel(sessions));

    for session = 1:numel(sessions)
        referenceDV{session} = 0;
        for area = 1:numel(referenceAreas)
            fileSession = sessions(session);
            if ~ismember(fileSession, areaSessions{strcmp(areas, referenceAreas{area})})
                fileSession = 11;   % session 1 saved as 11
            end
            R = load(fullfile(decoderFolder, sprintf('%s_%d_%s%s', ...
                subjects{subject}, fileSession, referenceAreas{area}, fileEnding)), 'log_odds_all');
            referenceDV{session} = referenceDV{session} + mean(R.log_odds_all(:,endWindow),2) / numel(referenceAreas);
        end
    end

    for area = 1:nAreas

        sessionVarianceSlope = nan(numel(sessions),1);

        % Pooled z-scored signed evidence and DV slopes of correct trials,
        % left (1) and right (2) choices
        pooledEvidence = {[], []};
        pooledSlope    = {[], []};

        residualCorrect = [];   % residual DV sign == choice, all trials

        for session = 1:numel(sessions)

            % S02 session 6: decoder output of 7PL, p47r and 8C not
            % reconstructed (MEG data missing for part of the recordings),
            % left out
            if strcmp(subjects{subject}, 'S02') && sessions(session) == 6 && ismember(areas{area}, {'7PL','p47r','8C'})
                continue
            end

            %% -----------------------------------------------------------
            % LOAD THIS SESSION
            %
            % log_odds_all     : trials x time, DV (> 0 = right choice)
            % Y_all            : -1 = left choice, +1 = right choice
            % stimulus_all     : -1 = left > right, +1 = right > left
            % evidence_strength: |mean right - mean left contrast|
            %% -----------------------------------------------------------

            fileSession = sessions(session);
            if ~ismember(fileSession, areaSessions{area})
                fileSession = 11;   % session 1 saved as 11
            end
            D = load(fullfile(decoderFolder, sprintf('%s_%d_%s%s', ...
                subjects{subject}, fileSession, areas{area}, fileEnding)));

            DV       = double(D.log_odds_all(:,1:200));
            choice   = sign(double(D.Y_all(:)));
            stimulus = double(D.stimulus_all(:));
            evidence = double(D.evidence_strength(:));


            %% -----------------------------------------------------------
            % CRITERION 1: SLOPE OF THE ACROSS-TRIAL VARIANCE (0.5 - 1 s)
            %% -----------------------------------------------------------

            choiceSlope = nan(1,2);

            for choiceSide = [-1 1]

                trials = choice == choiceSide;

                if sum(trials) >= 6
                    varianceTime = var(DV(trials,varianceWindow),0,1);
                    fitParameters = polyfit(timeSec(varianceWindow), varianceTime, 1);
                    choiceSlope((choiceSide+3)/2) = fitParameters(1);
                end

            end

            sessionVarianceSlope(session) = nanmean(choiceSlope);


            %% -----------------------------------------------------------
            % CRITERION 2: DV SLOPE VS EVIDENCE (correct trials)
            %% -----------------------------------------------------------

            trialSlope = nan(size(DV,1),1);

            for trial = 1:size(DV,1)
                fitParameters = polyfit(timeSec(slopeWindow)*1000, DV(trial,slopeWindow), 1);
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


            %% -----------------------------------------------------------
            % CRITERION 3: RESIDUAL DV (0.9 - 1 s) PREDICTS CHOICE
            %
            % PMd/M1 DV projected out within the session (regression
            % without intercept)
            %% -----------------------------------------------------------

            endDV      = mean(DV(:,endWindow),2);
            reference  = referenceDV{session};
            residualDV = endDV - (reference' * endDV) / (reference' * reference) * reference;

            predicted = sign(residualDV);
            residualCorrect = [residualCorrect; predicted(predicted ~= 0) == choice(predicted ~= 0)];

        end

        % Subject values
        varianceSlope(area,subject) = nanmean(sessionVarianceSlope);

        sideCorrelation = nan(1,2);
        for side = 1:2
            if numel(pooledEvidence{side}) >= 10
                sideCorrelation(side) = corr(pooledEvidence{side}, pooledSlope{side});
            end
        end
        evidenceScaling(area,subject) = nanmean(sideCorrelation);

        residualAccuracy(area,subject) = mean(residualCorrect);

    end

end


%% =======================================================================
% TABLE S1
%
% Mean and 95% confidence interval across participants
% (residual accuracy in %; reference areas without residual accuracy)
%% =======================================================================

measures     = {residualAccuracy * 100, varianceSlope, evidenceScaling};
measureNames = {'ResidualAccuracy', 'VarianceSlope', 'EvidenceScaling'};

tableS1 = table(tableAreas', 'VariableNames', {'ROI'});

for measure = 1:3

    meanValue = nan(numel(tableAreas),1);
    lowerCI   = nan(numel(tableAreas),1);
    upperCI   = nan(numel(tableAreas),1);

    for row = 1:numel(tableAreas)

        if measure == 1 && ismember(tableAreas{row}, referenceAreas)
            continue
        end

        values = measures{measure}(strcmp(areas, tableAreas{row}),:);
        values = values(isfinite(values));

        meanValue(row) = mean(values);
        halfWidth      = tinv(0.975, numel(values)-1) * std(values) / sqrt(numel(values));
        lowerCI(row)   = meanValue(row) - halfWidth;
        upperCI(row)   = meanValue(row) + halfWidth;

    end

    tableS1.([measureNames{measure} '_mean'])    = meanValue;
    tableS1.([measureNames{measure} '_CIlower']) = lowerCI;
    tableS1.([measureNames{measure} '_CIupper']) = upperCI;

end