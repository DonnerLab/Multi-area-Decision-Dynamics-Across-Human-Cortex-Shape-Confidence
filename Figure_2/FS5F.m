%% =======================================================================
% Figure
%% =======================================================================

clear; close all;

% Paths
scriptFolder  = fileparts(mfilename('fullpath'));
projectFolder = fileparts(fileparts(scriptFolder));

behaviorFile  = fullfile(projectFolder, 'data', 'behavior', 'data.mat');
modelFolder   = fullfile(projectFolder, 'data', 'model', 'results');
functionsPath = fullfile(projectFolder, 'functions');

addpath(functionsPath);

%% =======================================================================
% Parameters
%% =======================================================================

nSubjects    = 20;
nFolds       = 10;
nSimulations = 1;
modelNumber  = 4;

ffiValues = [0.05, 0.5, 0.8];
nFFI      = numel(ffiValues);

lag = 30;
timeWindow = 21:100;

All_sub_crosscorr = nan(nSubjects, 2 * lag + 1, nFFI);

%% =======================================================================
% Main loop
%% =======================================================================

load(behaviorFile);

for s = 1:nSubjects

    %% -------------------------------------------------------------------
    % Load behavioral data
    %% -------------------------------------------------------------------

    this_subject = [];

    for i = 1:3
        this_subject = [this_subject; cell2mat(All_Data(s,i))];
    end

    % Remove NaN trials
    this_subject(isnan(this_subject(:,1)), :) = [];

    data.t     = 10;
    data.left  = this_subject(:,2:11)';
    data.right = this_subject(:,12:21)';
    RT         = this_subject(:,24);

    %% -------------------------------------------------------------------
    % Load fitted model parameters
    %% -------------------------------------------------------------------

    modelFile = fullfile(modelFolder, ...
        ['Parameters_model_9_sbj_' num2str(s) '.mat']);

    load(modelFile);

    pars = Parameters(s,:,2);

    %% -------------------------------------------------------------------
    % Test different FFI values
    %% -------------------------------------------------------------------

    for ffi_n = 1:nFFI

        pars(5) = ffiValues(ffi_n);

        %% ---------------------------------------------------------------
        % Simulate model
        %% ---------------------------------------------------------------

        [Final_Confidence, Choices, ~, output, ~, ...
            Confidence_left, Confidence_right, ...
            DV_left, DV_right] = ...
            leaky_integrator_simulations(pars, data, modelNumber, nSimulations);

        DV_left_all  = DV_left;
        DV_right_all = DV_right;

        All_ch_model = Choices(:);
        All_ch_model(All_ch_model == 0) = -1;

        % RT proxy
        RT_model = abs(Final_Confidence(:));
        RT_model = 1 - RT_model;
        RT_model = RT_model + ...
            (mean(RT,'omitnan') - mean(RT_model,'omitnan'));

        %% ---------------------------------------------------------------
        % Train decoder from model DVs
        %% ---------------------------------------------------------------

        id_cv = cvpartition(size(All_ch_model,1), 'KFold', nFolds);
        Model_decoder = nan(size(All_ch_model,1), size(DV_left_all,1), 2);

        for side_n = 1:2
            for t = 1:size(DV_left_all,1)

                if side_n == 1
                    X_all = DV_right_all(t,:)';
                    Y_all = All_ch_model;
                else
                    X_all = DV_left_all(t,:)';
                    Y_all = -All_ch_model;
                end

                Prob_all = [];
                model_ch_dec = [];

                for cv = 1:nFolds

                    %% Training / test split
                    id_train = training(id_cv, cv);
                    id_test  = test(id_cv, cv);

                    X = X_all(id_train,:);
                    Y = Y_all(id_train);

                    %% Balance classes
                    id_left  = find(Y == -1);
                    id_right = find(Y == 1);

                    if numel(id_left) > numel(id_right)
                        keep_left  = id_left(randperm(numel(id_left), numel(id_right)));
                        id_all = [keep_left; id_right];
                    else
                        keep_right = id_right(randperm(numel(id_right), numel(id_left)));
                        id_all = [id_left; keep_right];
                    end

                    X = X(id_all,:);
                    Y = Y(id_all);

                    %% Normalize training data
                    [~, mu, sigma] = zscore(X);

                    %% Logistic decoder
                    Mdl = fitclinear( ...
                        zscore(X), Y, ...
                        'Learner', 'logistic');

                    %% Predict on test set
                    X_test = (X_all(id_test,:) - mu) ./ sigma;

                    Prob_fit = (X_test' .* Mdl.Beta) + Mdl.Bias;

                    Prob_all = [Prob_all; Prob_fit'];
                    model_ch_dec = [model_ch_dec; Choices(id_test)'];

                end

                Model_decoder(:,t,side_n) = Prob_all;

            end
        end

        %% ---------------------------------------------------------------
        % Winning vs losing DV
        %% ---------------------------------------------------------------

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

        %% ---------------------------------------------------------------
        % Cross-correlogram
        %% ---------------------------------------------------------------

        cross_corr = nan(size(P_left,1), 2 * lag + 1);

        for tr = 1:size(P_left,1)

            cross_corr(tr,:) = xcorr( ...
                zscore(DV_win(tr,timeWindow))', ...
                zscore(DV_lose(tr,timeWindow))', ...
                lag, ...
                'normalized');

        end

        All_sub_crosscorr(s,:,ffi_n) = mean(cross_corr,1,'omitnan');

    end
end

%% =======================================================================
% Plot
%% =======================================================================

figure;
colors = autumn(nFFI);

for i = 1:nFFI

    All_corrgram = All_sub_crosscorr(:,:,i);

    shadedErrorBar( ...
        -300:10:300, ...
        mean(All_corrgram,1,'omitnan'), ...
        std(All_corrgram,0,1,'omitnan') ./ sqrt(nSubjects), ...
        'lineProps', {'-', 'color', colors(i,:)} );

    hold on

end

xlabel('Time lag (ms)');
ylabel('Correlation');
set(gca,'fontsize',15);