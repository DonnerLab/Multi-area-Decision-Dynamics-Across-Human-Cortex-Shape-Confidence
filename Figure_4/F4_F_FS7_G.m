%% ============================================================
% PMd-PROJECTED NUMBER OF WINNING AREAS
%  Sliding Window Regression (Confidence + RT)
%% ============================================================

clear; close all; clc;
rng(42);

%% ============================================================
% PATHS
%% ============================================================

scriptFolder  = fileparts(mfilename('fullpath'));
projectFolder = fileparts(fileparts(scriptFolder));

dataPath = fullfile(projectFolder, 'data', 'decoder_output', 'all_areas');

%% ============================================================
% SETTINGS
%% ============================================================

subjects = {'S01','S04','S06','S07','S12','S15','S17','S18', ...
            'S11','S05','S19','S20','S22','S10','S02','S08', ...
            'S24','S21','S23','S25'};

nPerm = 10000;

%% ============================================================
% TIME SETTINGS
%% ============================================================

winLength = 35; %350 ms
stepSize  = 10; %100 ms steps

winStarts = 1:stepSize:(200 - winLength + 1);
nTime     = numel(winStarts);


timeSec = ((winStarts + winLength/2 - 1) * 10 - 500) / 1000;

%% ============================================================
% AREA SETS
%% ============================================================

% top 25 areas
set2 = {'AIP','FEF','5mv','23c','i6-8','24dd','6mp','PFt','IP2', ...
        '8Av','p24pr','7PC','23d','PF','6ma','55b','7AL','24dv', ...
        'DVT','PFm','31a','RSC','5L','V6','8C'};

% Pmd/M1 to project out
refAreas = {'6a','6d','4'};

areaSets = {set2(~ismember(set2, refAreas)),};

setNames = {'top 25 ROIs'};

nSets = 1;

%% ============================================================
% STORAGE
%% ============================================================

betasConfTime = cell(nSets, nTime);
betasRTTime   = cell(nSets, nTime);


%% ============================================================
% MAIN LOOP
%% ============================================================

for setIdx = 1

    areas = areaSets{setIdx};

    for subjIdx = 1:numel(subjects)

        subj = subjects{subjIdx};
        fprintf('Processing %s | Set %d\n', subj, setIdx);

        subjData   = struct();
        validAreas = {};
        validRef   = {};

        %% ----------------------------------------------------
        % LOAD ALL AREA FILES
        %% ----------------------------------------------------

        allAreas = unique([areas refAreas]);

        for areaIdx = 1:numel(allAreas)

            area = allAreas{areaIdx};

            files = dir(fullfile(dataPath, ...
                sprintf('%s_*_%s_CHOICEdecoder_alltrials_all.mat', subj, area)));

            if isempty(files)
                continue
            end

            logOdds = [];
            choice  = [];
            conf    = [];
            rt      = [];
            evid    = [];
            sess    = [];

            for f = 1:numel(files)

                D = load(fullfile(files(f).folder, files(f).name));
                X = D.log_odds_all(:,1:200);

                logOdds = [logOdds; X];
                choice  = [choice; D.Y_all(:)];
                conf    = [conf; D.confidence_all(:) - 1];
                rt      = [rt; D.rts_all(:)];
                evid    = [evid; D.evidence_strength(:)];
                sess    = [sess; f * ones(size(X,1),1)];
            end

            areaName = matlab.lang.makeValidName(area);

            subjData.(areaName).logodds = logOdds;
            subjData.(areaName).choice  = choice;
            subjData.(areaName).conf    = conf;
            subjData.(areaName).rt      = rt;
            subjData.(areaName).evid    = evid;
            subjData.(areaName).sess    = sess;

            if ismember(area, areas)
                validAreas{end+1} = areaName;
            end

            if ismember(area, refAreas)
                validRef{end+1} = areaName;
            end
        end

        if isempty(validAreas) || isempty(validRef)
            continue
        end

        %% ----------------------------------------------------
        % COMMON SESSIONS ACROSS ALL AREAS
        %% ----------------------------------------------------

        allNames = [validAreas validRef];

        sessAll = cellfun(@(x) unique(subjData.(x).sess), ...
            allNames, 'UniformOutput', false);

        commonSess = sessAll{1};

        for i = 2:numel(sessAll)
            commonSess = intersect(commonSess, sessAll{i});
        end

        if isempty(commonSess)
            continue
        end

        %% ----------------------------------------------------
        % RESTRICT TO COMMON SESSIONS
        %% ----------------------------------------------------

        for i = 1:numel(allNames)

            A = allNames{i};
            keep = ismember(subjData.(A).sess, commonSess);
            fn = fieldnames(subjData.(A));

            for j = 1:numel(fn)
                subjData.(A).(fn{j}) = subjData.(A).(fn{j})(keep,:);
            end
        end

        ref = validRef{1};

        choice     = subjData.(ref).choice;
        conf       = subjData.(ref).conf;
        rt         = subjData.(ref).rt;
        evid       = subjData.(ref).evid;
        sess       = subjData.(ref).sess;
        choiceSign = sign(choice);

        %% ====================================================
        % SLIDING WINDOW ANALYSIS
        %% ====================================================

        for t = 1:nTime

            idxT = winStarts(t):(winStarts(t) + winLength - 1);
            
            % calculate number of winning areas in this time window
            
            numWin = compute_NumWin_timepointwise( ...
                subjData, validAreas, validRef, idxT, choiceSign, sess);

            % zscore per session
            [numWinZ, evidZ, rtZ] = zscore_by_session(numWin, evid, rt, sess);

            keep = rt > 0;

            %% Confidence regression model
            mdlConf = fitglm( ...
                table(evidZ(keep), numWinZ(keep), conf(keep)), ...
                'Var3 ~ Var1 + Var2', ...
                'Distribution', 'binomial');

            %% RT regression model
            mdlRT = fitlm( ...
                table(evidZ(keep), numWinZ(keep), rtZ(keep)), ...
                'Var3 ~ Var1 + Var2');

            betasConfTime{setIdx, t} = [ ...
                betasConfTime{setIdx, t}; ...
                mdlConf.Coefficients.Estimate(2:3)'];

            betasRTTime{setIdx, t} = [ ...
                betasRTTime{setIdx, t}; ...
                mdlRT.Coefficients.Estimate(2:3)'];
        end

    
    end
end



%% ============================================================
% SLIDING WINDOW TIMECOURSE  CONFIDENCE
%% ============================================================

figure('Color','w','Position',[100 100 1200 300]);

for setIdx = 1

    X = nan(nTime, numel(subjects));

    for t = 1:nTime
        B = betasConfTime{setIdx, t};

        if ~isempty(B)
            X(t,1:size(B,1)) = B(:,2);
        end
    end

    sigMask = cluster_perm_1samp(X, nPerm, 0.05);

    meanB = mean(X, 2, 'omitnan');
    semB  = std(X, [], 2, 'omitnan') ./ sqrt(sum(~isnan(X),2));

    subplot(1, nSets, setIdx)
    hold on

    shadedErrorBar(timeSec, meanB, semB, 'lineprops', 'k');

    ySig = max(meanB) * 0.9;

    scatter(timeSec(sigMask), ...
        repmat(ySig, sum(sigMask), 1), ...
        30, 'r', 'filled');

    yline(0, 'k--');
    xline(0, 'k--');
    xline(1, 'k--');

    xlabel('Time from stimulus onset (s)');
    ylabel('\beta Number Winning Areas');
    title(setNames{setIdx});

    set(gca, 'FontSize', 15);
    box off
end

sgtitle(' Confidence´regression ');

%% ============================================================
% SLIDING WINDOW TIMECOURSE  REACTION TIME
%% ============================================================

figure('Color','w','Position',[100 200 1200 300]);

for setIdx = 1

    X = nan(nTime, numel(subjects));

    for t = 1:nTime
        B = betasRTTime{setIdx, t};

        if ~isempty(B)
            X(t,1:size(B,1)) = B(:,2);
        end
    end

    sigMask = cluster_perm_1samp(X, nPerm, 0.05);

    meanB = mean(X, 2, 'omitnan');
    semB  = std(X, [], 2, 'omitnan') ./ sqrt(sum(~isnan(X),2));

    subplot(1, nSets, setIdx)
    hold on

    shadedErrorBar(timeSec, meanB, semB, 'lineprops', 'k');

    ySig = max(abs(meanB)) * 1.5;

    scatter(timeSec(sigMask), ...
        repmat(ySig, sum(sigMask), 1), ...
        30, 'r', 'filled');

    yline(0, 'k--');
    xline(0, 'k--');
    xline(1, 'k--');

    xlabel('Time from stimulus onset (s)');
    ylabel('\beta Number Winning Areas');
    title(setNames{setIdx});

    set(gca, 'FontSize', 15);
    box off
end

sgtitle(' Reaction Time regression');