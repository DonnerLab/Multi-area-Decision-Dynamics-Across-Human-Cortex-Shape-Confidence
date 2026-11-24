%% =======================================================================
% Figure S6A-B
%
% Across-trial correlation between the average contrast (sample mean) of
% the left and right stimulus, after sorting trials by stimulus category
% (right > left, left > right) and choice
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

% Path to model results (behavior of each subject)
modelFolder = fullfile(projectFolder, 'data', 'model', 'results');

% Path to functions
addpath(fullfile(projectFolder, 'functions'));


%% =======================================================================
% SETTINGS
%% =======================================================================

nSubjects     = 20;
nPermutations = 10000;   % permutations of the permutation tests

% Participant shown in the scatter plots (same as Figure S1C)
exampleSubject = 1;


% Colors of the stimulus categories: right > left, left > right
categoryColors = [0.98 0.80 0.30; 0.90 0.40 0.10];

% subjects x category (1 = right > left, 2 = left > right) x choice
% (1 = left, 2 = right)
correlation_left_right = nan(nSubjects,2,2);


%% =======================================================================
% MAIN LOOP ACROSS SUBJECTS
%% =======================================================================

for subject = 1:nSubjects

    load(fullfile( ...
        modelFolder, ...
        ['All_data_model_9_sbj_' num2str(subject) '.mat']));

    % Sample mean of each side, choice of the participant (1 = right, 0 = left)
    meanLeft   = mean(All_left,2);
    meanRight  = mean(All_right,2);
    realChoice = Real_data(:,1);

    rightGreater = meanRight > meanLeft;

    for choice = 1:2

        choiceTrials = realChoice == (choice - 1);

        for category = 1:2

            if category == 1
                trials = choiceTrials & rightGreater;
            else
                trials = choiceTrials & ~rightGreater;
            end

            correlation_left_right(subject,category,choice) = ...
                corr(meanLeft(trials), meanRight(trials));

        end


        %% ===============================================================
        % SCATTER PLOT OF THE EXAMPLE PARTICIPANT (left part of S6A / S6B)
        %% ===============================================================

        if subject == exampleSubject

            figure(choice)
            subplot(1,2,1)
            hold on

            % Background: stimulus categories
            fill([0.2 0.8 0.2], [0.2 0.8 0.8], categoryColors(1,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
            fill([0.2 0.8 0.8], [0.2 0.2 0.8], categoryColors(2,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none');

            scatter(meanLeft(choiceTrials), meanRight(choiceTrials), 15, 'k', 'filled', ...
                'MarkerFaceAlpha', 0.3);

            % Linear fit within each stimulus category
            for category = 1:2

                if category == 1
                    trials = choiceTrials & rightGreater;
                else
                    trials = choiceTrials & ~rightGreater;
                end

                fitParameters = polyfit(meanLeft(trials), meanRight(trials), 1);
                fitX = linspace(0.2, 0.8, 100);

                plot(fitX, polyval(fitParameters, fitX), 'color', categoryColors(category,:), 'linewidth', 4);

            end

            title(sprintf('pooled r = %.2f', corr(meanLeft(choiceTrials), meanRight(choiceTrials))), ...
                'FontWeight', 'normal')
            axis square
            xlim([0.2 0.8])
            ylim([0.2 0.8])
            xlabel('Single-trial sample mean left')
            ylabel('Single-trial sample mean right')
            set(gca,'fontsize',15)

        end

    end

end


%% =======================================================================
% CORRELATION COEFFICIENTS OF ALL PARTICIPANTS (right part of S6A / S6B)
%% =======================================================================

choiceNames = {'Left choice', 'Right choice'};

for choice = 1:2

    figure(choice)
    subplot(1,2,2)
    hold on

    pValues = nan(1,2);

    for category = 1:2

        pValues(category) = permutationTest( ...
            correlation_left_right(:,category,choice), ...
            zeros(nSubjects,1), ...
            nPermutations);

        bar(category, mean(correlation_left_right(:,category,choice)), ...
            'FaceColor', categoryColors(category,:));

        scatter(category + (rand(nSubjects,1)*2 - 1)*0.05, correlation_left_right(:,category,choice), ...
            80, 'k', 'filled', 'MarkerFaceAlpha', 0.5);

    end

    xticks([1 2])
    xticklabels({'Right > Left','Left > Right'})
    ylabel('Correlation coefficient')
    title({choiceNames{choice}; sprintf('p = %.2g,  p = %.2g', pValues(1), pValues(2))}, ...
        'FontWeight', 'normal')
    set(gca,'fontsize',15)

end
