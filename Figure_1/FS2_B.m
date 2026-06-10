%% =======================================================================
% Model comparison: Mean Squared Error (MSE)
%% =======================================================================

clearvars

%% Paths

% Folder where this script is located
scriptFolder = fileparts(mfilename('fullpath'));

% Go up two folders
repoRoot = (fileparts(scriptFolder));

% Standard project paths
dataFolder      = fullfile(repoRoot, 'data');
modelFolder     = fullfile(dataFolder, 'model', 'results');
functionsFolder = fullfile(repoRoot, 'functions');

% Add functions (including VBA toolbox if stored here)
addpath(genpath(functionsFolder));

%% Settings

nModels   = 19;
nSubjects = 20;

colors_model = autumn(10);
colors_model = repmat(colors_model, 2, 1);

%% Compute MSE

figure(3)
hold on

for model_n = 1:nModels

    All_MSE = nan(nSubjects,1);

    for s = 1:nSubjects

        load(fullfile(modelFolder, ...
            sprintf('All_data_model_%d_sbj_%d.mat', model_n, s)));

        %% Extract variables

        Side = double(mean(All_right,2) > mean(All_left,2));
        Side(Side == 0) = -1;

        Conf_data_sbj = Real_data(:,2)';
        All_GM_sbj    = Real_data(:,4);

        All_left_sbj  = All_left;
        All_right_sbj = All_right;

        % Predicted probabilities
        All_prob_sbj = abs(All_Conf);

        % Model choice probabilities
        All_ch_model = sum(All_prob_sbj(:,1:2),2);

        % Data choices and RT
        All_RT      = Real_data(:,3); 
        All_ch_data = Real_data(:,1);

        %% Contrast fluctuations around generative mean

        All_right = All_right_sbj - (0.5 + All_GM_sbj ./ 2);
        All_left  = All_left_sbj  - (0.5 - All_GM_sbj ./ 2);

        %% Confidence reports

        All_Conf_Data = Conf_data_sbj' - 1;

        % Simulated choice/confidence reports from model
        nTrials = size(All_ch_model,1);
        All_cc = zeros(nTrials,100,4);

        for trial = 1:nTrials
            All_cc(trial,:,:) = mnrnd(1, All_prob_sbj(trial,:), 100);
        end

        %% Trial categories: model

        id_1_1 = find(All_cc(:,:,1) == 1);
        id_1_0 = find(All_cc(:,:,2) == 1);
        id_0_1 = find(All_cc(:,:,3) == 1);
        id_0_0 = find(All_cc(:,:,4) == 1);

        %% Trial categories: data

        id_1_1_d = find(All_ch_data == 1 & All_Conf_Data == 1);
        id_1_0_d = find(All_ch_data == 1 & All_Conf_Data == 0);
        id_0_1_d = find(All_ch_data == 0 & All_Conf_Data == 1);
        id_0_0_d = find(All_ch_data == 0 & All_Conf_Data == 0);

        %% Confidence kernels

        All_right_all = All_right;
        All_left_all  = All_left;

        datapoints_1 = [];
        datapoints_2 = [];
        datapoints_3 = [];
        datapoints_4 = [];

        model_1 = [];
        model_2 = [];
        model_3 = [];
        model_4 = [];

        for samp = 1:10

            All_right = repmat(All_right_all(:,samp),100,4);
            All_left  = repmat(All_left_all(:,samp),100,4);

            % Data

            datapoints_1 = [datapoints_1; ...
                mean([mean(mean(All_right(id_1_1_d),1),3); ...
                      mean(mean(All_left(id_0_1_d),1),3)])];

            datapoints_2 = [datapoints_2; ...
                mean([mean(mean(All_right(id_0_1_d),1),3); ...
                      mean(mean(All_left(id_1_1_d),1),3)])];

            datapoints_3 = [datapoints_3; ...
                mean([mean(mean(All_right(id_1_0_d),1),3); ...
                      mean(mean(All_left(id_0_0_d),1),3)])];

            datapoints_4 = [datapoints_4; ...
                mean([mean(mean(All_right(id_0_0_d),1),3); ...
                      mean(mean(All_left(id_1_0_d),1),3)])];

            % Model

            model_1 = [model_1; ...
                mean([mean(mean(All_right(id_1_1),1),3); ...
                      mean(mean(All_left(id_0_1),1),3)])];

            model_2 = [model_2; ...
                mean([mean(mean(All_right(id_0_1),1),3); ...
                      mean(mean(All_left(id_1_1),1),3)])];

            model_3 = [model_3; ...
                mean([mean(mean(All_right(id_1_0),1),3); ...
                      mean(mean(All_left(id_0_0),1),3)])];

            model_4 = [model_4; ...
                mean([mean(mean(All_right(id_0_0),1),3); ...
                      mean(mean(All_left(id_1_0),1),3)])];

        end

        %% Mean squared error

        MSE_1 = mean((datapoints_1 - model_1).^2);
        MSE_2 = mean((datapoints_2 - model_2).^2);
        MSE_3 = mean((datapoints_3 - model_3).^2);
        MSE_4 = mean((datapoints_4 - model_4).^2);

        All_MSE(s) = MSE_1 + MSE_2 + MSE_3 + MSE_4;

    end

    %% Plot

    bar(model_n, mean(All_MSE), ...
        'FaceColor', colors_model(model_n,:), ...
        'FaceAlpha', 1);

    errorbar(model_n, ...
        mean(All_MSE), ...
        std(All_MSE) / sqrt(nSubjects), ...
        'Color', colors_model(model_n,:), ...
        'LineStyle', 'none');

end

ylabel('MSE')

set(gca, ...
    'FontSize', 15, ...
    'Box', 'off', ...
    'TickDir', 'out');