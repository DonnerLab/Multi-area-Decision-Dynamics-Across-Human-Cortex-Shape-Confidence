%% =======================================================================
% Figure 1C and Figure S1E-F
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

% Path to model results (behavior + model predictions of each subject)
modelFolder = fullfile(projectFolder, 'data', 'model', 'results');

% Path to functions
addpath(fullfile(projectFolder, 'functions'));


%% =======================================================================
% SETTINGS
%% =======================================================================

nSubjects     = 20;
nSamples      = 10;
nBins         = 7;        % delta-contrast bins of the psychometric curve
nPermutations = 10000;   % permutations of the permutation tests

% Participants with a primacy index <= this value are shown in red
% (participants without a primacy effect)
primacyThreshold = 0.05;



P_rightChoice       = nan(nSubjects,nBins);
P_rightChoice_model = nan(nSubjects,nBins);
binCenters          = nan(nBins,1);

AUC_all = nan(nSubjects,nSamples);


%% =======================================================================
% MAIN LOOP ACROSS SUBJECTS
%% =======================================================================

for subject = 1:nSubjects

    %% ===================================================================
    % LOAD SUBJECT DATA
    %
    % All_right / All_left : trials x samples contrasts
    % All_ch_model         : model probability of a right choice
    % Real_data            : [choice, confidence, RT, generative mean]
    %% ===================================================================

    load(fullfile( ...
        modelFolder, ...
        ['All_data_model_9_sbj_' num2str(subject) '.mat']));

    % Mean contrast difference (right - left) and sample-wise difference
    evidence           = mean(All_right,2) - mean(All_left,2);
    contrastDifference = All_right - All_left;

    % Participant choice (1 = right, 0 = left) and model choice
    realChoice  = Real_data(:,1);
    modelChoice = All_ch_model(:,1);

    % Remove trials without a response
    noResponse = isnan(evidence) | isnan(realChoice);

    evidence(noResponse)             = [];
    contrastDifference(noResponse,:) = [];
    realChoice(noResponse)           = [];
    modelChoice(noResponse)          = [];


    %% ===================================================================
    % PSYCHOMETRIC CURVE (Figure 1C)
    %
    %% ===================================================================

    quantiles = round(quantile(abs(evidence),4),2);

    [~, edges, evidenceBin] = histcounts( ...
        evidence, ...
        [flip(-quantiles) quantiles]);

    for bin = 1:nBins

        P_rightChoice(subject,bin)       = nanmean(realChoice(evidenceBin == bin));
        P_rightChoice_model(subject,bin) = nanmean(modelChoice(evidenceBin == bin));
        binCenters(bin)                  = mean(edges(bin:bin+1));

    end


    %% ===================================================================
    % PSYCHOPHYSICAL KERNEL (Figure S1E)
    %
      %% ===================================================================

    rightGreater = evidence > 0;
    leftGreater  = evidence < 0;

    AUC_rightGreater = nan(1,nSamples);
    AUC_leftGreater  = nan(1,nSamples);

    for sample = 1:nSamples

        [~,~,~,AUC_rightGreater(sample)] = perfcurve( ...
            realChoice(rightGreater), ...
            contrastDifference(rightGreater,sample), ...
            1);

        [~,~,~,AUC_leftGreater(sample)] = perfcurve( ...
            realChoice(leftGreater), ...
            contrastDifference(leftGreater,sample), ...
            1);

    end

    AUC_all(subject,:) = mean([AUC_leftGreater; AUC_rightGreater],1) - 0.5;

end


%% =======================================================================
% PSYCHOPHYSICAL KERNEL STATISTICS
%
%% =======================================================================

pValues_kernel = nan(1,nSamples);

for sample = 1:nSamples
    pValues_kernel(sample) = permutationTest( ...
        AUC_all(:,sample), ...
        zeros(nSubjects,1), ...
        nPermutations);
end

significantSamples = double(fdr_bh(pValues_kernel,0.05));
significantSamples(significantSamples == 0) = nan;


%% =======================================================================
% PRIMACY INDEX (Figure S1F)
%

%% =======================================================================

firstHalf  = sum(AUC_all(:,1:5),2);
secondHalf = sum(AUC_all(:,6:10),2);

primacyIndex = (firstHalf - secondHalf) ./ (abs(firstHalf) + abs(secondHalf));

pValue_primacy = permutationTest( ...
    primacyIndex, ...
    zeros(nSubjects,1), ...
    nPermutations);

% Participants without a primacy effect (shown in red)
noPrimacy = primacyIndex <= primacyThreshold;


%% =======================================================================
% FIGURE 1 (Figure 1C) - PSYCHOMETRIC CURVE
%% =======================================================================

figure(1)
hold on

% Data: group average +/- SEM
e = errorbar( ...
    binCenters, ...
    nanmean(P_rightChoice,1), ...
    nanstd(P_rightChoice,0,1) ./ sqrt(nSubjects), ...
    'color', 'k', ...
    'linewidth', 3, ...
    'linestyle', 'none');
e.CapSize = 0;

scatter( ...
    binCenters, ...
    nanmean(P_rightChoice,1), ...
    100, ...
    'k', ...
    'filled');

% Model: group average +/- SEM
shadedErrorBar( ...
    binCenters, ...
    nanmean(P_rightChoice_model,1), ...
    nanstd(P_rightChoice_model,0,1) ./ sqrt(nSubjects), ...
    'lineProps', {'-','color','k','linewidth',3});

box off
xticks(binCenters)
ylabel('Probability right choice')
xlabel('\Delta Contrast (R-L)')
set(gca,'fontsize',15)


%% =======================================================================
% FIGURE 2 (Figure S1E) - PSYCHOPHYSICAL KERNEL
%% =======================================================================

figure(2)
hold on

% Single participants (no primacy effect in red, plotted on top)
for subject = [find(~noPrimacy)' find(noPrimacy)']

    if noPrimacy(subject)
        lineColor = [0.85 0.10 0.10];
    else
        lineColor = [0.5 0.5 0.5];
    end

    plot(1:nSamples, AUC_all(subject,:), 'color', lineColor, 'linewidth', 2);

end

% Group average
plot(1:nSamples, nanmean(AUC_all,1), 'color', 'k', 'linewidth', 5);

% Significant samples (p < 0.05, FDR-corrected)
scatter( ...
    1:nSamples, ...
    significantSamples .* (max(AUC_all(:)) * 1.1), ...
    80, ...
    'k', ...
    'square', ...
    'filled');

box off
xticks([1 nSamples])
xlim([0 nSamples+1])
ylabel('ROC index - 0.5')
xlabel('Sample number')
set(gca,'fontsize',15)


%% =======================================================================
% FIGURE 3 (Figure S1F) - PRIMACY INDEX
%% =======================================================================

figure(3)
hold on

bar(1, mean(primacyIndex), 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'k');

for subject = 1:nSubjects

    if noPrimacy(subject)
        dotColor = [0.85 0.10 0.10];
    else
        dotColor = [0 0 0];
    end

    scatter( ...
        1 + (rand(1) - 0.5) * 0.08, ...
        primacyIndex(subject), ...
        80, ...
        dotColor, ...
        'filled');

end

yline(0,'--k');

xlim([0.5 1.5])
xticks([])
ylabel('Primacy index')
title(sprintf('p = %.2g', pValue_primacy))
box off
set(gca,'FontSize',15)
