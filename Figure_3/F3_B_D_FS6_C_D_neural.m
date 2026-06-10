%% ============================================================
% PMd/M1 HEMISPHERE DV REGRESSION
% Sliding Timecourse + Asymmetry Analysis
%% ============================================================

clear;
close all;
clc;
rng(42);

%% ============================================================
% PATHS
%% ============================================================

scriptFolder  = fileparts(mfilename('fullpath'));
projectFolder = fileparts(fileparts(scriptFolder));

dataPath = fullfile(projectFolder, ...
    'data', 'decoder_output', '6a6d4');

%% ============================================================
% SETTINGS
%% ============================================================

subjects = { ...
    'S01','S04','S06','S07','S12','S15','S17','S18', ...
    'S11','S05','S19','S20','S22','S10','S02','S08', ...
    'S24','S21','S23','S25'};

nSubjects = numel(subjects);
nSessions = 6;
nPerm     = 10000;

%% ============================================================
% TIME WINDOWS
%% ============================================================

winLength = 35;   % ~350 ms
stepSize  = 10;   % ~100 ms

winStarts = 1:stepSize:(201 - winLength + 1);
nTime     = numel(winStarts);

% Decoder timeline: -500 ms to +1500 ms


timeSec = ((winStarts + winLength/2 - 1) * 10 - 500) ./ 1000;

% stimulus offset window 
finalWindow = 134:168;

%% ============================================================
% STORAGE
%% ============================================================

betasConfTime = cell(1, nTime);
betasRTTime   = cell(1, nTime);

% Confidence:
% [evidence beta, asymmetry (abs(win)-abs(lose)), p-value evidence]
betasConfLast = [];

% RT:
% [evidence beta, asymmetry(abs(win)-abs(lose)),]
betasRTLast = [];

%% ============================================================
% MAIN LOOP
%% ============================================================

for s = 1:nSubjects

    subject = subjects{s};
    fprintf('Processing %s\n', subject);

    DV_left  = [];
    DV_right = [];

    choice   = [];
    conf     = [];
    rt       = [];
    evidence = [];
    sessID   = [];

    %% --------------------------------------------------------
    % LOAD DATA
    %% --------------------------------------------------------

    for ses = 1:nSessions

        for hemi = 1:2

            if hemi == 1
                hemiLabel = 'l';
            else
                hemiLabel = 'r';
            end

            fileName = fullfile(dataPath, ...
                sprintf('%s_SES%d_%s_CHOICEdecoder_6A6D4_all.mat', ...
                subject, ses, hemiLabel));

            if ~isfile(fileName)
                continue
            end

            D = load(fileName);
            % decoder log odds
            X = D.log_odds_all(:,1:201);

            if hemi == 1
                DV_left = [DV_left; X];

                % behavioral variables
                choice   = [choice; D.Y_all(:)];
                conf     = [conf; D.confidence_all(:) - 1];
                rt       = [rt; D.rts_all(:)];
                evidence = [evidence; D.evidence_strength(:)];
                sessID   = [sessID; ses .* ones(size(X,1),1)];
            else
                DV_right = [DV_right; X];
            end
        end
    end

    if isempty(DV_left) || isempty(DV_right)
        continue
    end

    choiceSign = sign(choice);

    %% ========================================================
    % SLIDING WINDOW REGRESSION
    %% ========================================================

    for t = 1:nTime

        idx = winStarts(t):(winStarts(t) + winLength - 1);

        % for this time window define winning and losing DV
        [DV_win, DV_lose] = compute_winlose_DV(...
            DV_left, DV_right, idx, choiceSign);

        [DV_win, DV_lose, evidence_z, rt_z] = ...
            zscore_by_session(...
            DV_win, DV_lose, evidence, rt, sessID);

        keep = rt > 0;

        %% ---------------- CONFIDENCE regression ----------------

        mdl = fitglm(...
            table(...
            evidence_z(keep), ...
            DV_win(keep), ...
            DV_lose(keep), ...
            conf(keep)), ...
            'Var4 ~ Var1 + Var2 + Var3', ...
            'Distribution', 'binomial');

        betasConfTime{t} = [ ...
            betasConfTime{t}; ...
            mdl.Coefficients.Estimate(2:4)'];

        %% ---------------- RT regression----------------

        mdl = fitlm(...
            table(...
            evidence_z(keep), ...
            DV_win(keep), ...
            DV_lose(keep), ...
            rt_z(keep)), ...
            'Var4 ~ Var1 + Var2 + Var3');

        betasRTTime{t} = [ ...
            betasRTTime{t}; ...
            mdl.Coefficients.Estimate(2:4)'];
    end

    %% ========================================================
    % stimulus offset ANALYSIS
    %% ========================================================

    [DV_win, DV_lose] = compute_winlose_DV(...
        DV_left, DV_right, finalWindow, choiceSign);

    [DV_win, DV_lose, evidence_z, rt_z] = ...
        zscore_by_session(...
        DV_win, DV_lose, evidence, rt, sessID);

    keep = rt >= 0;

    %% ---------------- CONFIDENCE ----------------

    mdl = fitglm(...
        table(...
        evidence_z(keep), ...
        DV_win(keep), ...
        DV_lose(keep), ...
        conf(keep)), ...
        'Var4 ~ Var1 + Var2 + Var3', ...
        'Distribution', 'binomial');

    b = mdl.Coefficients.Estimate(2:4);
    p = mdl.Coefficients.pValue(2:4);

    asymmetry = abs(b(2)) - abs(b(3));

    betasConfLast = [ ...
        betasConfLast; ...
        [b(1), asymmetry, p(1)]];

    %% ---------------- RT ----------------

    mdl = fitlm(...
        table(...
        evidence_z(keep), ...
        DV_win(keep), ...
        DV_lose(keep), ...
        rt_z(keep)), ...
        'Var4 ~ Var1 + Var2 + Var3');

    b = mdl.Coefficients.Estimate(2:4);
    asymmetry = abs(b(2)) - abs(b(3));

    betasRTLast = [ ...
        betasRTLast; ...
        [b(1), asymmetry]];
end

%% ============================================================
% FILTER SUBJECTS
%% ============================================================

% Keep subjects with:
% evidence beta > 0
% AND significant evidence effect

evidenceBeta = betasConfLast(:,1);
pEvidence    = betasConfLast(:,3);

keepSubjects = evidenceBeta > 0 & pEvidence < 0.05;

fprintf('Keeping %d / %d subjects\n', ...
    sum(keepSubjects), numel(keepSubjects));

betasConfLast = betasConfLast(keepSubjects, 1:2);

%% ============================================================
% PLOTS
%% ============================================================

plot_three_beta_timecourse(...
    betasConfTime, timeSec, nPerm, 'Confidence');

plot_three_beta_timecourse(...
    betasRTTime, timeSec, nPerm, 'Reaction Time');

plot_last_window(...
    betasConfLast, 'Confidence (final window)');

plot_last_window(...
    betasRTLast, 'Reaction Time (final window)');

%% ============================================================
% FUNCTIONS
%% ============================================================

function [DVwin, DVlose] = compute_winlose_DV(...
    DV_left, DV_right, idx, choiceSign)

leftDV  = mean(DV_left(:, idx), 2);
rightDV = mean(DV_right(:, idx), 2);

DVwin  = nan(size(leftDV));
DVlose = nan(size(leftDV));

rightChoice = choiceSign == 1;

DVwin(rightChoice)  = leftDV(rightChoice);
DVlose(rightChoice) = rightDV(rightChoice);

DVwin(~rightChoice)  = rightDV(~rightChoice);
DVlose(~rightChoice) = leftDV(~rightChoice);

end

function [a, b, c, d] = zscore_by_session(...
    a, b, c, d, sessID)

allSessions = unique(sessID);

for s = 1:numel(allSessions)

    idx = sessID == allSessions(s);

    if sum(idx) > 2
        a(idx) = zscore(a(idx));
        b(idx) = zscore(b(idx));
        c(idx) = zscore(c(idx));
        d(idx) = zscore(d(idx));
    end
end

end

function plot_three_beta_timecourse(...
    betaCell, timeSec, nPerm, plotTitle)

nTime = numel(betaCell);
nSub  = max(cellfun(@(x) size(x,1), betaCell));

X = nan(nTime, 3, nSub);

for t = 1:nTime

    B = betaCell{t};

    if ~isempty(B)
        X(t, 1:size(B,2), 1:size(B,1)) = permute(B, [2 1]);
    end
end

figure('Color','w','Position',[100 100 1000 300]);
hold on

colors = [ 0 0 0;...
    0.00 0.45 0.90; ...
    ];

for r = 2:3

    M = squeeze(X(:, r, :));

    mu = mean(M, 2, 'omitnan');
    se = std(M, [], 2, 'omitnan') ./ ...
        sqrt(sum(~isnan(M), 2));

    shadedErrorBar(...
        timeSec, mu, se, ...
        'lineProps', {'-','Color', colors(r-1,:), 'LineWidth', 2});

    sigMask = cluster_perm_1samp(M, nPerm, 0.05);

    ySig = max(abs(mean(M, 2, 'omitnan'))) * 1.5 + 0.02*(r-2);

    scatter(...
        timeSec(sigMask), ...
        repmat(ySig, sum(sigMask), 1), ...
        25, colors(r-1,:), 'filled');
end

box off
set(gca, 'FontSize', 14)
xlabel('Time (s)')
ylabel('Beta weight')
title(plotTitle)

end

function plot_last_window(B, plotTitle)

figure('Color','w','Position',[400 400 500 300]);

labels = { ...
    'Evidence strength', ...
    'Asymmetry |win| - |lose|'};

for k = 1:2

    subplot(1,2,k)

    y = B(:,k);

    bar(1, mean(y), 0.5, ...
        'FaceColor', 'k', ...
        'FaceAlpha', 0.2);
    hold on

    scatter(...
        1 + (rand(size(y))*0.1 - 0.05), ...
        y, ...
        100, 'k', 'filled', ...
        'MarkerFaceAlpha', 0.5);

    yline(0, 'k--')

    p = permutationTest(...
        y, zeros(size(y)), 10000);

    text(...
        1, mean(y), ...
        sprintf('p = %.3g', p), ...
        'HorizontalAlignment', 'center');

    title(labels{k})
    xticks([])
    box off
end

sgtitle(plotTitle)

end