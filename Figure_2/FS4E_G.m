%% =======================================================================
% Figure S4E-G
%
%
% S4E: decoder accuracy over time, choice and confidence decoder
% S4F: cosine similarity between the weight vectors of the two decoders
% S4G: weights assigned to each PC by the two decoders at stimulus
%      offset (example subject, session and cross-validation fold)
%
% Each subject: average across cross-validation folds, then sessions
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

% Decoder output
decoderFolder = fullfile(projectFolder, 'data', 'decoder_output', 'all_areas');
fileSuffix    = '_CHOICE_CONFIDENCE_decoder_all_vx__win_beta_net_kfold.mat';

% Path to functions
addpath(fullfile(projectFolder, 'functions'));


%% =======================================================================
% SETTINGS
%% =======================================================================

subjects = { ...
    'S01','S04','S06','S07','S12','S15','S17','S18','S11','S05', ...
    'S19','S20','S22','S10','S02','S08','S24','S21','S23','S25'};

nSubjects = numel(subjects);

% Example of Figure S4G
exampleSubject = 'S01';
exampleSession = 11;
exampleFold    = 1;

choiceColor     = [0.00 0.45 0.74];
confidenceColor = [0.85 0.33 0.10];


%% =======================================================================
% LOAD ACCURACY AND COSINE SIMILARITY (subjects x time)
%% =======================================================================

for subject = 1:nSubjects

    files = dir(fullfile(decoderFolder, [subjects{subject} '_SES*' fileSuffix]));

    sessionAccuracy_choice     = [];
    sessionAccuracy_confidence = [];
    sessionCosine              = [];

    for file = 1:numel(files)

        % accuracy_choice, accuracy_conf, cos_sim: folds x time
        D = load(fullfile(decoderFolder, files(file).name), ...
            'accuracy_choice', 'accuracy_conf', 'cos_sim');

        sessionAccuracy_choice     = [sessionAccuracy_choice;     nanmean(double(D.accuracy_choice),1)];
        sessionAccuracy_confidence = [sessionAccuracy_confidence; nanmean(double(D.accuracy_conf),1)];
        sessionCosine              = [sessionCosine;              nanmean(double(D.cos_sim),1)];

    end

    accuracy_choice(subject,:)     = nanmean(sessionAccuracy_choice,1);
    accuracy_confidence(subject,:) = nanmean(sessionAccuracy_confidence,1);
    cosineSimilarity(subject,:)    = nanmean(sessionCosine,1);

end

% Time axis: -0.5 to 1.5 s from stimulus onset
timeSec = linspace(-0.5, 1.5, size(accuracy_choice,2));


%% =======================================================================
% FIGURE 1 (Figure S4E) - DECODER ACCURACY
%% =======================================================================

figure(1)
hold on

shadedErrorBar(timeSec, nanmean(accuracy_choice,1), ...
    nanstd(accuracy_choice,0,1) ./ sqrt(nSubjects), ...
    'lineProps', {'-','color',choiceColor,'linewidth',2});

shadedErrorBar(timeSec, nanmean(accuracy_confidence,1), ...
    nanstd(accuracy_confidence,0,1) ./ sqrt(nSubjects), ...
    'lineProps', {'-','color',confidenceColor,'linewidth',2});

yline(0.5, '--k');
xline(0, '-', 'Color', [0.5 0.5 0.5]);
xline(1, '-', 'Color', [0.5 0.5 0.5]);

xlim([-0.5 1.5])
xticks([-0.5 0 1 1.5])
xlabel('Time from stimulus onset (s)')
ylabel({'Decoder accuracy'; '(fraction correct)'})
title('Choice (blue), confidence (orange)', 'FontWeight', 'normal')
box off
set(gca,'FontSize',15)


%% =======================================================================
% FIGURE 2 (Figure S4F) - COSINE SIMILARITY OF THE DECODER WEIGHTS
%% =======================================================================

figure(2)
hold on

shadedErrorBar(timeSec, nanmean(cosineSimilarity,1), ...
    nanstd(cosineSimilarity,0,1) ./ sqrt(nSubjects), ...
    'lineProps', {'-','color',[0.2 0.6 0.3],'linewidth',2});

yline(0, '--k');
xline(0, '-', 'Color', [0.5 0.5 0.5]);
xline(1, '-', 'Color', [0.5 0.5 0.5]);

xlim([-0.5 1.5])
xticks([-0.5 0 1 1.5])
xlabel('Time from stimulus onset (s)')
ylabel('Cosine similarity')
box off
set(gca,'FontSize',15)


%% =======================================================================
% FIGURE 3 (Figure S4G) - DECODER WEIGHTS PER PC AT STIMULUS OFFSET
%
% beta_choice_pc / beta_conf_pc: PCs x time x folds (NaN after the last
% PC used); n_pcs_all: folds x time, number of PCs used
%% =======================================================================

example = load(fullfile(decoderFolder, [exampleSubject '_SES' num2str(exampleSession) fileSuffix]), ...
    'beta_choice_pc', 'beta_conf_pc', 'n_pcs_all');

% Time point closest to stimulus offset (1 s)
[~, offsetTime] = min(abs(timeSec - 1));

nPCs = double(example.n_pcs_all(exampleFold,offsetTime));

weights_choice     = double(example.beta_choice_pc(1:nPCs,offsetTime,exampleFold));
weights_confidence = double(example.beta_conf_pc(1:nPCs,offsetTime,exampleFold));

maxWeight = max(abs([weights_choice; weights_confidence]));

figure(3)

subplot(2,1,1)
bar(1:nPCs, weights_choice, 'FaceColor', choiceColor, 'EdgeColor', 'none');
xlim([0.5 nPCs+0.5])
ylim([-1.1 1.1] * maxWeight)
ylabel('Weight')
title('Choice', 'FontWeight', 'normal')
box off
set(gca,'FontSize',15)

subplot(2,1,2)
bar(1:nPCs, weights_confidence, 'FaceColor', confidenceColor, 'EdgeColor', 'none');
xlim([0.5 nPCs+0.5])
ylim([-1.1 1.1] * maxWeight)
xlabel('PCs')
ylabel('Weight')
title('Confidence', 'FontWeight', 'normal')
box off
set(gca,'FontSize',15)
