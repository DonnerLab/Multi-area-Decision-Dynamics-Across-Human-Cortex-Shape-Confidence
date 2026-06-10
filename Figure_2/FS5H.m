%% =======================================================================
% Figure S5H correaltion
%% =======================================================================

clear; close all;

% Paths
scriptFolder  = fileparts(mfilename('fullpath'));
projectFolder = fileparts(fileparts(scriptFolder));

dataFolder = fullfile(projectFolder, 'data', 'decoder_output', 'bilateral', 'stimulus');
paramFile  = fullfile(projectFolder, 'data', 'decoder_output', 'Parameters_all.mat');

%% =======================================================================
% Subjects
%% =======================================================================

SubjectID = {'S01','S02','S04','S05','S07','S08','S11','S12', ...
             'S15','S17','S18','S10','S06','S19','S20','S21','S22','S23','S24','S25'};

nSub = numel(SubjectID);
nSes = 6;

Mean_slope = nan(nSub, nSes);

%% =======================================================================
% Loop over subjects and sessions
%% =======================================================================

for s = 1:nSub
    for ses = 1:nSes

        baseName = sprintf('%s_choice_both_session%d', SubjectID{s}, ses);

        pcaFile     = fullfile(dataFolder, [baseName '_PCA_PMd_M1_100Hz.mat']);
        decoderFile = fullfile(dataFolder, [baseName '_PMd_M1_decoder_output_100hz.mat']);

        if ~isfile(pcaFile) || ~isfile(decoderFile)
            continue
        end

        % Load data
        load(pcaFile)      % contains Choices
        load(decoderFile)  % contains Prob_Performance

        %% ----------------------------------------------------------------
        % Decoder performance (logistic transform)
        %% ----------------------------------------------------------------
        Perf = 1 ./ (1 + exp(Prob_Performance));
        Perf = permute(Perf, [3,2,1]);   % trials × time × ?

        %% ----------------------------------------------------------------
        % Compute slope per trial
        %% ----------------------------------------------------------------
        nTrials = size(Perf,1);
        slopes  = nan(nTrials,1);

        timeVec = (0.2:0.01:1)';  % matches indices 71:151

        for tr = 1:nTrials
            y = Perf(tr,71:151)';
            p = polyfit(timeVec, y, 1);
            slopes(tr) = abs(p(1));
        end

        Mean_slope(s,ses) = mean(slopes,'omitnan');

    end
end

%% =======================================================================
% Load model parameters
%% =======================================================================

load(paramFile)   % contains All_pars
Leak = All_pars(:,1);

%% =======================================================================
% Correlation plot
%% =======================================================================

meanSlope = mean(Mean_slope,2,'omitnan');

figure;
scatter(meanSlope, Leak, 100, 'k', 'filled', 'MarkerFaceAlpha', 0.5);
hold on;

[r,p] = corr(meanSlope, Leak, 'Type','spearman');

% linear fit
fitLine = polyfit(meanSlope, Leak', 1);
xFit = linspace(min(meanSlope), max(meanSlope), 100);
yFit = polyval(fitLine, xFit);
plot(xFit, yFit, 'k', 'LineWidth', 1.5);

title(sprintf('Spearman r = %.3f, p = %.3f', r, p));
xlabel('Neural DV slope');
ylabel('Leak');
set(gca,'fontsize',15);