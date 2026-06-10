%% ========================================================================
% Figure 2F, 2H and S4E - model DV
% 
% 
% Load behavioral + model data and simulate model decoder correlations
%% ========================================================================

clear all
close all

%% ========================================================================
% Paths
%% ========================================================================

% Folder where this script is located
scriptFolder = fileparts(mfilename('fullpath'));

% Go up two folders
projectFolder = fileparts(fileparts(scriptFolder));

% Main folders
behaviorFolder  = fullfile(projectFolder, 'data', 'behavior');
modelFolder     = fullfile(projectFolder, 'data', 'model', 'results');
functionsFolder = fullfile(projectFolder, 'functions','model_fit');

% Add required paths
addpath(functionsFolder)

%% ========================================================================
% Initialization
%% ========================================================================

nSubjects = 20;
nFolds    = 10;
nTime     = 100;

All_sub_corr  = nan(nSubjects,1);
All_corr_RT   = nan(nSubjects,1);
All_corr_conf = nan(nSubjects,1);

%% ========================================================================
% Loop across subjects
%% ========================================================================

for s = 1:nSubjects

    %% --------------------------------------------------------------------
    % Load behavioral data
    %% --------------------------------------------------------------------

    load(fullfile(behaviorFolder, 'data.mat'))

    this_subject = [];

    for i = 1:3
        this_subject = [this_subject; cell2mat(All_Data(s,i))];
    end

    % Remove NaN trials
    id_nan = isnan(this_subject(:,1));
    this_subject(id_nan,:) = [];

    data.t     = 10;
    data.left  = this_subject(:,2:11)';
    data.right = this_subject(:,12:21)';
    RT         = this_subject(:,24);

    %% --------------------------------------------------------------------
    % Load fitted model parameters
    %% --------------------------------------------------------------------

    load(fullfile(modelFolder, ...
        ['Parameters_model_9_sbj_' num2str(s) '.mat']))

    pars = Parameters(s,:,2);

    %% --------------------------------------------------------------------
    % Simulate best-fitting model (Model 4)
    %% --------------------------------------------------------------------

    [Final_Confidence, Choices, Confidence, output, bound, ...
        Confidence_left, Confidence_right, DV_left, DV_right] = ...
        leaky_integrator_simulations(pars, data, 4, 1);

    %% --------------------------------------------------------------------
    % Extract model variables
    %% --------------------------------------------------------------------

    DV_left_all  = DV_left;
    DV_right_all = DV_right;

    % Stimulus evidence
    Ev = mean(data.right - data.left, 1);

    % Stimulus side
    Side = double(sign(Ev));
    Side(Side == -1) = 0;
    Side = Side';

    % Model choices
    All_ch_model = Choices(:);

    % Confidence reports
    All_conf_model = output(:);

    % RT proxy from model confidence
    RT_model = abs(Final_Confidence(:));
    RT_model = 1 - RT_model;
    RT_model = RT_model + (nanmean(RT) - nanmean(RT_model));

    %% --------------------------------------------------------------------
    % Train decoder on model DVs
    %% --------------------------------------------------------------------

    id_cv = cvpartition(size(All_ch_model,1), 'KFold', nFolds);

    Model_decoder = [];

    for side_n = 1:2
        for t = 1:nTime

            model_ch_dec   = [];
            model_conf_dec = [];
            model_RT_dec   = [];

            %% Select decoder input

            if side_n == 1
                X_all = DV_right_all(t,:)';
            else
                X_all = DV_left_all(t,:)';
            end

            %% Choice labels

            All_ch_model_dec = double(All_ch_model);
            All_ch_model_dec(All_ch_model_dec == 0) = -1;

            if side_n == 1
                Y_all = All_ch_model_dec;
            else
                Y_all = -All_ch_model_dec;
            end

            Y_all(Y_all == 0) = -1;

            %% Cross-validation

            Prob_all = [];

            for cv = 1:nFolds

                % Train/test split
                id_train = training(id_cv, cv);
                id_test  = test(id_cv, cv);

                X = X_all(id_train,:);
                Y = Y_all(id_train);

                %% Balance classes

                id_left  = find(Y == -1);
                id_right = find(Y == 1);

                if length(id_left) > length(id_right)

                    id_all = [ ...
                        id_right;
                        id_left(randperm(length(id_left), length(id_right)))
                        ];

                else

                    id_all = [ ...
                        id_left;
                        id_right(randperm(length(id_right), length(id_left)))
                        ];
                end

                X = X(id_all,:);
                Y = Y(id_all);

                %% Normalize training data

                [~, mu, sigma] = zscore(X);

                %% Train logistic decoder

                Mdl = fitclinear( ...
                    zscore(X), ...
                    Y, ...
                    'Learner', 'logistic');

                %% Predict test data

                norm_x_test = (X_all(id_test,:) - mu) ./ sigma;

                Y_fit = predict(Mdl, norm_x_test);

                % Log odds
                Prob_fit = ...
                    (norm_x_test' .* Mdl.Beta) + Mdl.Bias;

                %% Concatenate outputs

                Prob_all = [Prob_all; Prob_fit'];

                model_ch_dec   = [model_ch_dec; All_ch_model(id_test)];
                model_conf_dec = [model_conf_dec; All_conf_model(id_test)];
                model_RT_dec   = [model_RT_dec; RT_model(id_test)];

            end

            Model_decoder(:,t,side_n) = Prob_all;

        end
    end

    %% ====================================================================
    % Winning vs losing DV assignment
    %% ====================================================================

    P_right = Model_decoder(:,:,1);
    P_left  = Model_decoder(:,:,2);

    DV_win  = nan(size(P_left));
    DV_lose = nan(size(P_left));

    for tr = 1:size(P_left,1)

        if model_ch_dec(tr) == 0
            DV_win(tr,:)  = P_left(tr,:);
            DV_lose(tr,:) = P_right(tr,:);
        else
            DV_win(tr,:)  = P_right(tr,:);
            DV_lose(tr,:) = P_left(tr,:);
        end
    end

    %% ====================================================================
    % Cross-correlation between winning and losing DV
    %% ====================================================================

    cross_corr = [];

    for tr = 1:size(P_left,1)

        lag = 10;

        par = xcorr( ...
            zscore(DV_win(tr,21:100))', ...
            zscore(DV_lose(tr,21:100))', ...
            lag, ...
            'normalized');

        cross_corr(tr,:) = par;
    end

    corr_all = mean(cross_corr,2);

    %% ====================================================================
    % Subject-level summary statistics
    %% ====================================================================

    All_sub_corr(s) = mean(corr_all);

    All_corr_RT(s) = corr( ...
        model_RT_dec, ...
        corr_all, ...
        'type', 'pearson');

    Conf_fit = glmfit( ...
        zscore(corr_all), ...
        model_conf_dec, ...
        'binomial', ...
        'link', 'logit');

    All_corr_conf(s) = Conf_fit(2);

end

%% ========================================================================
% Plot: mean cross-correlation
%% ========================================================================

figure(1)

bar(1, mean(All_sub_corr), ...
    'facecolor', 'k', ...
    'facealpha', 0.5)

hold on

for s = 1:nSubjects
    scatter( ...
        1 + (rand(1)*2 - 1) * 0.05, ...
        All_sub_corr(s), ...
        100, ...
        'k', ...
        'filled', ...
        'markerfacealpha', 0.5)
end

ylabel('Correlation coefficient')
set(gca, 'fontsize', 15)

%% ========================================================================
% Plot: correlation with model confidence
%% ========================================================================

figure(2)

bar(1, mean(All_corr_conf), ...
    'facecolor', 'k', ...
    'facealpha', 0.5)

hold on

for s = 1:nSubjects
    scatter( ...
        1 + (rand(1)*2 - 1) * 0.05, ...
        All_corr_conf(s), ...
        100, ...
        'k', ...
        'filled', ...
        'markerfacealpha', 0.5)
end

ylabel('Correlation with model confidence')
set(gca, 'fontsize', 15)

%% ========================================================================
% Plot: correlation with model RT
%% ========================================================================

figure(3)

bar(1, mean(All_corr_RT), ...
    'facecolor', 'k', ...
    'facealpha', 0.5)

hold on

for s = 1:nSubjects
    scatter( ...
        1 + (rand(1)*2 - 1) * 0.05, ...
        All_corr_RT(s), ...
        100, ...
        'k', ...
        'filled', ...
        'markerfacealpha', 0.5)
end

ylabel('Correlation with model RT')
set(gca, 'fontsize', 15)