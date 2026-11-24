%% =======================================================================
% Figure S2C
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

% Path to model results
modelFolder = fullfile(projectFolder, 'data', 'model', 'results');


%% =======================================================================
% SETTINGS
%% =======================================================================

nModels      = 19;
nSubjects    = 20;
nSimulations = 1000;   

modelColors = repmat(autumn(10), 2, 1);

meanSquaredError = nan(nSubjects,nModels);


%% =======================================================================
% MEAN SQUARED ERROR OF EACH MODEL AND SUBJECT
%% =======================================================================

for model = 1:nModels
    for subject = 1:nSubjects

        load(fullfile(modelFolder, ...
            ['All_data_model_' num2str(model) '_sbj_' num2str(subject) '.mat']));

        % Participant choice (1 = right, 0 = left) and confidence (1 = high, 0 = low)
        realChoice     = Real_data(:,1);
        realConfidence = Real_data(:,2) - 1;

        % Contrast fluctuations around the generative mean of each side
        generativeMean   = Real_data(:,4);
        fluctuationRight = All_right - (0.5 + generativeMean ./ 2);
        fluctuationLeft  = All_left  - (0.5 - generativeMean ./ 2);

        % Report of each trial (trials x 4):
        % 1 = right/high, 2 = right/low, 3 = left/high, 4 = left/low
        % Data: observed report; model: number of simulations (out of
        % nSimulations) with each report
        report_data = double([ ...
            realChoice == 1 & realConfidence == 1, ...
            realChoice == 1 & realConfidence == 0, ...
            realChoice == 0 & realConfidence == 1, ...
            realChoice == 0 & realConfidence == 0]);

        report_model = mnrnd(nSimulations, abs(All_Conf) ./ sum(abs(All_Conf),2));

        % Average contrast fluctuation of each report (4 reports x 10 samples)
        kernelRight_data  = (report_data'  * fluctuationRight) ./ sum(report_data,1)';
        kernelLeft_data   = (report_data'  * fluctuationLeft)  ./ sum(report_data,1)';
        kernelRight_model = (report_model' * fluctuationRight) ./ sum(report_model,1)';
        kernelLeft_model  = (report_model' * fluctuationLeft)  ./ sum(report_model,1)';

        % The 4 confidence kernels (right and left choices averaged):
        % selected/high, non-selected/high, selected/low, non-selected/low
        kernels_data = [ ...
            mean([kernelRight_data(1,:); kernelLeft_data(3,:)],1); ...
            mean([kernelRight_data(3,:); kernelLeft_data(1,:)],1); ...
            mean([kernelRight_data(2,:); kernelLeft_data(4,:)],1); ...
            mean([kernelRight_data(4,:); kernelLeft_data(2,:)],1)];

        kernels_model = [ ...
            mean([kernelRight_model(1,:); kernelLeft_model(3,:)],1); ...
            mean([kernelRight_model(3,:); kernelLeft_model(1,:)],1); ...
            mean([kernelRight_model(2,:); kernelLeft_model(4,:)],1); ...
            mean([kernelRight_model(4,:); kernelLeft_model(2,:)],1)];

        % Mean across samples, summed across the 4 kernels
        meanSquaredError(subject,model) = sum(mean((kernels_data - kernels_model).^2, 2));

    end
end


%% =======================================================================
% FIGURE (Figure S2C)
%% =======================================================================

figure
hold on

for model = 1:nModels

    bar(model, nanmean(meanSquaredError(:,model)), ...
        'FaceColor', modelColors(model,:));

    errorbar(model, ...
        nanmean(meanSquaredError(:,model)), ...
        nanstd(meanSquaredError(:,model)) / sqrt(nSubjects), ...
        'Color', modelColors(model,:), ...
        'LineStyle', 'none');

end

ylabel('Mean squared error')
xlabel('Model')
set(gca, 'XTick', 1:nModels, 'FontSize', 15, 'Box', 'off', 'TickDir', 'out');
