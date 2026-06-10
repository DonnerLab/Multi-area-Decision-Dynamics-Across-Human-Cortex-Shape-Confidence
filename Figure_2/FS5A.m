%% =======================================================================
% Model-based variance and mean DV analysis (Figure S4A)
%% =======================================================================

clear
close all

%% =======================================================================
% Paths
%% =======================================================================

% Folder where this script is located
scriptFolder = fileparts(mfilename('fullpath'));

% Go up two folders
projectFolder = fileparts(fileparts(scriptFolder));

% Data folders
behaviorFolder  = fullfile(projectFolder, 'data', 'behavior');
modelFolder     = fullfile(projectFolder, 'data', 'model', 'results');
functionsFolder = fullfile(projectFolder, 'Figure 5', 'functions');

addpath(functionsFolder)

%% =======================================================================
% Main loop across subjects
%% =======================================================================

for s = 1:20

    %% -------------------------------------------------------------------
    % Load behavioral data
    %% -------------------------------------------------------------------

    load(fullfile(behaviorFolder, 'data.mat'))

    this_subject = [];

    for i = 1:3
        this_subject = [this_subject; cell2mat(All_Data(s, i))];
    end

    id_nan = isnan(this_subject(:,1));
    this_subject(id_nan, :) = [];

    data.t     = 10;
    data.left  = this_subject(:, 2:11)';
    data.right = this_subject(:,12:21)';
    RT         = this_subject(:,24);

    %% -------------------------------------------------------------------
    % Load best-fitting model parameters
    %% -------------------------------------------------------------------

    load(fullfile(modelFolder, ...
        ['Parameters_model_9_sbj_', num2str(s), '.mat']))

    pars = Parameters(s,:,2);

    %% -------------------------------------------------------------------
    % Simulate best model
    %% -------------------------------------------------------------------

    [Final_Confidence, Choices, Confidence, output, bound, ...
        Confidence_left, Confidence_right, ...
        DV_left, DV_right, P_cross] = ...
        leaky_integrator_simulations(pars, data, 4, 1);

    %% -------------------------------------------------------------------
    % Basic variables
    %% -------------------------------------------------------------------

    Evidence = mean(data.right - data.left, 1);

    % stimulus strength
    model_Ev = abs(Evidence(:));

    % stimulus side: left = 0, right = 1
    Side = sign(Evidence);
    Side(Side == -1) = 0;
    Side = Side';

    % model choices
    All_ch_model = Choices(:);

    %% -------------------------------------------------------------------
    % Stop integration after bound crossing
    %% -------------------------------------------------------------------

    for tr = 1:size(Confidence_left, 2)

        id_cross = (P_cross(:, tr) == 1);

        DV_left(id_cross, tr)  = Confidence_left(tr);
        DV_right(id_cross, tr) = Confidence_right(tr);

    end

    %% -------------------------------------------------------------------
    % Original model DV (pre-decoder)
    %% -------------------------------------------------------------------

    DV_win_orig  = nan(size(DV_left))';
    DV_lose_orig = nan(size(DV_left))';

    for tr = 1:size(DV_win_orig, 1)

        if All_ch_model(tr) == 1
            DV_win_orig(tr, :)  = DV_right(:, tr);
            DV_lose_orig(tr, :) = DV_left(:, tr);
        else
            DV_win_orig(tr, :)  = DV_left(:, tr);
            DV_lose_orig(tr, :) = DV_right(:, tr);
        end

    end

    %% -------------------------------------------------------------------
    % Model confidence and RT proxy
    %% -------------------------------------------------------------------

    All_conf_model = output(:);

    RT_model = abs(Final_Confidence(:));
    RT_model = 1 - RT_model;
    RT_model = RT_model + (nanmean(RT) - nanmean(RT_model));

    %% -------------------------------------------------------------------
    % Train linear decoder from model DV
    %% -------------------------------------------------------------------

    DV_left_all  = DV_left;
    DV_right_all = DV_right;

    id_cv = cvpartition(size(All_ch_model,1), 'KFold', 10);
    Model_decoder = [];

    for side_n = 1:2
        for t = 1:size(DV_left_all,1)

            if side_n == 1
                X_all = DV_right_all(t,:)';
            else
                X_all = DV_left_all(t,:)';
            end

            All_ch_model_dec = double(All_ch_model);
            All_ch_model_dec(All_ch_model_dec == 0) = -1;

            if side_n == 1
                Y_all = All_ch_model_dec;
            else
                Y_all = -All_ch_model_dec;
            end

            Y_all(Y_all == 0) = -1;

            Prob_all = [];
            model_ch_dec   = [];
            model_conf_dec = [];
            model_RT_dec   = [];

            %% 10-fold CV

            for cv = 1:10

                id_train = training(id_cv, cv);
                id_test  = test(id_cv, cv);

                X = X_all(id_train, :);
                Y = Y_all(id_train);

                %% Upsample minority class

                id_left  = find(Y == -1);
                id_right = find(Y == 1);

                if numel(id_left) > numel(id_right)

                    id_keep = [ ...
                        id_right;
                        id_left(randperm(numel(id_left), numel(id_right)))];

                else

                    id_keep = [ ...
                        id_left;
                        id_right(randperm(numel(id_right), numel(id_left)))];

                end

                X = X(id_keep, :);
                Y = Y(id_keep);

                %% Train decoder

                [~, mu, sigma] = zscore(X);

                Mdl = fitclinear( ...
                    zscore(X), ...
                    Y, ...
                    'Learner', 'logistic');

                %% Test decoder

                X_test = (X_all(id_test,:) - mu) ./ sigma;

                Prob_fit = (X_test' .* Mdl.Beta) + Mdl.Bias;

                Prob_all = [Prob_all; Prob_fit'];

                model_ch_dec   = [model_ch_dec; All_ch_model(id_test)];
                model_conf_dec = [model_conf_dec; All_conf_model(id_test)];
                model_RT_dec   = [model_RT_dec; RT_model(id_test)];

            end

            Model_decoder(:, t, side_n) = Prob_all;

        end
    end

    %% -------------------------------------------------------------------
    % Decoder-derived DV (log-odds)
    %% -------------------------------------------------------------------

    P_right = Model_decoder(:,:,1);
    P_left  = Model_decoder(:,:,2);

    DV_win  = nan(size(P_left));
    DV_lose = nan(size(P_left));

    for tr = 1:size(P_left,1)

        if model_ch_dec(tr) == 1
            DV_win(tr,:)  = P_right(tr,:);
            DV_lose(tr,:) = P_left(tr,:);
        else
            DV_win(tr,:)  = P_left(tr,:);
            DV_lose(tr,:) = P_right(tr,:);
        end

    end

    %% -------------------------------------------------------------------
    % Correct trials only
    %% -------------------------------------------------------------------

    id_correct = find(Side == All_ch_model);

    for time_p = 1:20

        idx = (1:5) + 5 * (time_p - 1);

        % across-trial variance
        Var_win(time_p) = nanvar( ...
            nanmean(DV_win_orig(id_correct, idx), 2));

        Var_lose(time_p) = nanvar( ...
            nanmean(DV_lose_orig(id_correct, idx), 2));

        % mean decoder DV
        Mean_win(time_p) = mean( ...
            mean(DV_win(id_correct, idx), 2));

        Mean_lose(time_p) = mean( ...
            mean(DV_lose(id_correct, idx), 2));

    end

    %% Store subject results

    All_within_contra_sep(:, s) = Var_win(:);
    All_within_ipsi_sep(:, s)   = Var_lose(:);

    All_mean_contra_sep(:, s) = Mean_win(:);
    All_mean_ipsi_sep(:, s)   = Mean_lose(:);

end

%% =======================================================================
% Plot variance
%% =======================================================================

figure(1)

subplot(2,1,1)

shadedErrorBar( ...
    1:size(All_within_contra_sep,1), ...
    nanmean(All_within_contra_sep, 2), ...
    nanstd(All_within_contra_sep, 0, 2) ./ sqrt(size(All_within_contra_sep,2)), ...
    'lineProps', {'-', 'color', 'k'})

hold on

shadedErrorBar( ...
    1:size(All_within_ipsi_sep,1), ...
    nanmean(All_within_ipsi_sep, 2), ...
    nanstd(All_within_ipsi_sep, 0, 2) ./ sqrt(size(All_within_ipsi_sep,2)), ...
    'lineProps', {'-', 'color', 'b'})

xticks([1 8 24 30])
xticklabels([-0.5 0 1 1.5])

xlabel('Time from stimulus onset')
ylabel('Variance')
set(gca, 'FontSize', 15)

%% =======================================================================
% Plot mean decoder output
%% =======================================================================

subplot(2,1,2)

shadedErrorBar( ...
    1:size(All_mean_contra_sep,1), ...
    nanmean(All_mean_contra_sep, 2), ...
    nanstd(All_mean_contra_sep, 0, 2) ./ sqrt(size(All_mean_contra_sep,2)), ...
    'lineProps', {'-', 'color', 'k'})

hold on

shadedErrorBar( ...
    1:size(All_mean_ipsi_sep,1), ...
    nanmean(All_mean_ipsi_sep, 2), ...
    nanstd(All_mean_ipsi_sep, 0, 2) ./ sqrt(size(All_mean_ipsi_sep,2)), ...
    'lineProps', {'-', 'color', 'b'})

xticks([1 8 24 30])
xticklabels([-0.5 0 1 1.5])

xlabel('Time from stimulus onset')
ylabel('Mean DV')
set(gca, 'FontSize', 15)