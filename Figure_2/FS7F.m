%% =======================================================================
% Figure S7F
%
% Correlation across participants between the leak parameter of the
% best-fitting model and the slope of the neural DV (bilateral PMd/M1
% choice decoder) during stimulus viewing
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

% Bilateral PMd/M1 decoder output (100 Hz) and model results
decoderFolder = fullfile(projectFolder, 'data', 'decoder_output', 'bilateral', 'stimulus');
modelFolder   = fullfile(projectFolder, 'data', 'model', 'results');


%% =======================================================================
% SETTINGS
%% =======================================================================

% Subjects in the order of the model fits
subjects = { ...
    'S01','S02','S04','S05','S07','S08','S11','S12','S15','S17', ...
    'S18','S10','S06','S19','S20','S21','S22','S23','S24','S25'};

nSubjects = numel(subjects);
nSessions = 6;

% Samples 71:151 = 0.2 - 1 s from stimulus onset (10 ms each)
slopeSamples = 71:151;
slopeTime    = (0.2:0.01:1)';

meanSlope = nan(nSubjects,nSessions);
leak      = nan(nSubjects,1);


%% =======================================================================
% MAIN LOOP ACROSS SUBJECTS
%% =======================================================================

for subject = 1:nSubjects

    for session = 1:nSessions

        decoderFile = fullfile(decoderFolder, sprintf('%s_choice_both_session%d_PMd_M1_decoder_output_100hz.mat', ...
            subjects{subject}, session));

        if ~isfile(decoderFile)
            continue
        end

        % Prob_Performance: 1 x time x trials, decoder output
        load(decoderFile)

        % Probability of the decoded choice, trials x time
        decoderProbability = permute(1 ./ (1 + exp(Prob_Performance)), [3 2 1]);

        % Absolute slope of each trial
        trialSlope = nan(size(decoderProbability,1),1);

        for trial = 1:size(decoderProbability,1)
            fitParameters = polyfit(slopeTime, decoderProbability(trial,slopeSamples)', 1);
            trialSlope(trial) = abs(fitParameters(1));
        end

        meanSlope(subject,session) = nanmean(trialSlope);

    end

    % Leak of the best-fitting model
    load(fullfile(modelFolder, ['Parameters_model_9_sbj_' num2str(subject) '.mat']));
    leak(subject) = Parameters(subject,1,2);

end

meanSlope = nanmean(meanSlope,2);


%% =======================================================================
% CORRELATION AND FIGURE (Figure S7F)
%% =======================================================================

[rho, pValue] = corr(meanSlope, leak, 'Type', 'Spearman');

fitParameters = polyfit(meanSlope, leak, 1);
fitX = linspace(min(meanSlope), max(meanSlope), 100);

figure
hold on

scatter(meanSlope, leak, 100, 'k', 'filled', 'MarkerFaceAlpha', 0.5);
plot(fitX, polyval(fitParameters, fitX), 'k', 'LineWidth', 1.5);

title(sprintf('Rho = %.2f, p = %.2g', rho, pValue), 'FontWeight', 'normal')
xlabel('Neural DV slope')
ylabel('Leak \lambda')
box off
set(gca, 'fontsize', 15)
