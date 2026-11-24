%% =======================================================================
% Figure S1C-D
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


%% =======================================================================
% SETTINGS
%% =======================================================================

nSubjects = 20;

% Participant shown in Figure S1C
exampleSubject = 1;


green  = [0.60 0.80 0.55];
orange = [0.98 0.80 0.55];

% subjects x category (1 = right > left, 2 = left > right)
correlation_left_right = nan(nSubjects,2);


%% =======================================================================
% MAIN LOOP ACROSS SUBJECTS
%% =======================================================================

for subject = 1:nSubjects

    load(fullfile( ...
        modelFolder, ...
        ['All_data_model_9_sbj_' num2str(subject) '.mat']));

    % Contrast fluctuations around the generative mean of each side
    generativeMean   = Real_data(:,4);
    fluctuationRight = All_right - (0.5 + generativeMean ./ 2);
    fluctuationLeft  = All_left  - (0.5 - generativeMean ./ 2);

    % Stimulus category
    rightGreater = mean(All_right,2) >  mean(All_left,2);
    leftGreater  = mean(All_right,2) <= mean(All_left,2);

    % Correlation across all samples of all trials of each category
    right_rightGreater = fluctuationRight(rightGreater,:);
    left_rightGreater  = fluctuationLeft(rightGreater,:);
    right_leftGreater  = fluctuationRight(leftGreater,:);
    left_leftGreater   = fluctuationLeft(leftGreater,:);

    correlation_left_right(subject,1) = corr(left_rightGreater(:), right_rightGreater(:));
    correlation_left_right(subject,2) = corr(left_leftGreater(:),  right_leftGreater(:));


    %% ===================================================================
    % FIGURE 1 (Figure S1C)
    %% ===================================================================

    if subject == exampleSubject

        figure(1)

        subplot(1,2,1)
        scatter(right_rightGreater(:), left_rightGreater(:), 10, green, 'filled');
        title({'Right > Left'; sprintf('Rho: %.3f', correlation_left_right(subject,1))})
        xlabel('R excess contrast')
        ylabel('L excess contrast')
        axis square
        set(gca,'fontsize',15)

        subplot(1,2,2)
        scatter(right_leftGreater(:), left_leftGreater(:), 10, orange, 'filled');
        title({'Left > Right'; sprintf('Rho: %.3f', correlation_left_right(subject,2))})
        xlabel('R excess contrast')
        ylabel('L excess contrast')
        axis square
        set(gca,'fontsize',15)

    end

end


%% =======================================================================
% FIGURE 2 (Figure S1D) - 
%% =======================================================================

figure(2)
hold on

bar(1:2, mean(correlation_left_right,1), 'facecolor', 'k', 'facealpha', 0.5);

errorbar( ...
    1:2, ...
    mean(correlation_left_right,1), ...
    std(correlation_left_right,0,1) ./ sqrt(nSubjects), ...
    'color', 'k', ...
    'linestyle', 'none');

for category = 1:2
    for subject = 1:nSubjects
        scatter( ...
            category + (rand(1)*2 - 1) * 0.05, ...
            correlation_left_right(subject,category), ...
            100, ...
            'k', ...
            'filled', ...
            'markerfacealpha', 0.5);
    end
end

ylim([-1 1])
xticks([1 2])
xticklabels({'Right > Left','Left > Right'})
ylabel('Correlation coefficients')
set(gca,'fontsize',15)
