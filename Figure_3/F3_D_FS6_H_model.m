%% Clean Model GLM Analysis
% Figure 5B-F and S13
% Confidence / RT regression using winning vs losing decision variables

clear; close all; clc;

%% ================================================================
%% PATHS
%% ================================================================

% Automatic project path generation
scriptFolder  = fileparts(mfilename('fullpath'));
projectFolder = fileparts(fileparts(scriptFolder));

behaviorFile  = fullfile(projectFolder, 'data', 'behavior', 'data.mat');
modelFolder   = fullfile(projectFolder, 'data', 'model', 'results');
functionsPath = fullfile(projectFolder, 'Figure 5', 'functions');

addpath(functionsPath);

%% ================================================================
%% MODEL DEFINITIONS
%% ================================================================

model_names = {'model 4', 'model 17', 'model 18'};
model_ids   = [4, 17, 18];
sim_ids     = [4, 5, 7];

n_subjects = 20;
n_simulations = 1;

%% ================================================================
%% LOAD BEHAVIORAL DATA
%% ================================================================

load(behaviorFile, 'All_Data')

%% ================================================================
%% RUN ANALYSIS
%% ================================================================

for model_idx = 2

    model_number = model_ids(model_idx);
    simulation_id = sim_ids(model_idx);

    fprintf('\nRunning %s\n', model_names{model_idx})

    %% ------------------------------------------------------------
    %% LOAD PARAMETERS ACROSS SUBJECTS
    %% ------------------------------------------------------------

    if model_number > 17
        all_pars = nan(n_subjects, 12);
    else
        all_pars = nan(n_subjects, 11);
    end

    pooled_subject_data = [];

    for s = 1:n_subjects

        %% Concatenate subject behavioral sessions
        subject_data = [];

        for block = 1:3
            subject_data = [subject_data; cell2mat(All_Data(s, block))];
        end

        subject_data(isnan(subject_data(:,1)), :) = [];

        pooled_subject_data = [pooled_subject_data; subject_data];

        %% Load fitted model parameters
        parameter_file = fullfile( ...
            modelFolder, ...
            sprintf('Parameters_model_%d_sbj_%d.mat', model_number, s));

        load(parameter_file, 'Parameters')
        all_pars(s,:) = Parameters(s,:,2);
    end

    %% Mean parameters across subjects
    pars = nanmean(all_pars, 1);

    %% Manual inhibition increase for these models
    if ismember(simulation_id, [5 6 7])
        pars(5) = 1;
    end

    %% ------------------------------------------------------------
    %% PREPARE INPUT DATA
    %% ------------------------------------------------------------

    data.t = 10;
    data.left  = pooled_subject_data(:,2:11)';
    data.right = pooled_subject_data(:,12:21)';
    RT = pooled_subject_data(:,24);

    %% ------------------------------------------------------------
    %% SIMULATE MODEL
    %% ------------------------------------------------------------

    [Final_Confidence, Choices, Confidence, output, bound, ...
        Confidence_left, Confidence_right] = ...
        leaky_integrator_simulations(pars, data, simulation_id, n_simulations);

    %% ------------------------------------------------------------
    %% EXTRACT MODEL VARIABLES
    %% ------------------------------------------------------------

    % Evidence strength
    evidence = abs(mean(data.right - data.left, 1));
    evidence = repmat(evidence, [1, 1, n_simulations]);
    evidence = evidence(:);

    % Model choices
    model_choices = Choices(:);

    % Final decision variables
    dv_left  = Confidence_left(:);
    dv_right = Confidence_right(:);

    % Confidence report
    model_confidence = output(:);

    % RT proxy
    model_rt = abs(Final_Confidence(:));
    model_rt = 1 - model_rt;
    model_rt = model_rt + (nanmean(RT) - nanmean(model_rt));

    %% ------------------------------------------------------------
    %% WINNING / LOSING DVs
    %% ------------------------------------------------------------

    dv_win  = nan(size(dv_left));
    dv_lose = nan(size(dv_left));

    for tr = 1:length(dv_left)

        if model_choices(tr) == 1
            dv_win(tr)  = dv_right(tr);
            dv_lose(tr) = dv_left(tr);
        else
            dv_win(tr)  = dv_left(tr);
            dv_lose(tr) = dv_right(tr);
        end
    end

    %% ------------------------------------------------------------
    %% RIDGE LOGISTIC REGRESSION FOR CONFIDENCE
    %% ------------------------------------------------------------

    X = [abs(evidence), dv_win, dv_lose];
    Xz = zscore(X);
    Ycat = categorical(model_confidence);

    mdl_conf = fitclinear( ...
        Xz, Ycat, ...
        'Learner', 'logistic', ...
        'Regularization', 'ridge', ...
        'Lambda', 1e-3);

    glm_weights_conf = mdl_conf.Beta;

    %% ------------------------------------------------------------
    %% RIDGE LINEAR REGRESSION FOR RT
    %% ------------------------------------------------------------

    glm_weights_rt = ridge(model_rt, X, 1e-3, 0);
    glm_weights_rt = glm_weights_rt(2:end);

    %% Store
    All_GLM_weights_Conf = glm_weights_conf';
    All_GLM_weights_RT   = glm_weights_rt';

    %% ============================================================
    %% PLOT RESULTS
    %% ============================================================

    for plot_type = 1

        if plot_type == 1
            all_glm = All_GLM_weights_Conf;
        else
            all_glm = All_GLM_weights_RT;
        end

        colors = [ ...
            0.0000    0.4470    0.7410
            235 67 52
            235 128 52
            0.4660    0.6740    0.1880] ./ [1 1 1; 256 256 256; 256 256 256; 1 1 1];

        figure(model_idx)
        clf
        sgtitle(model_names{model_idx})

        %% --------------------------------------------------------
        %% Panel 1: Winning vs Losing DV weights
        %% --------------------------------------------------------

        subplot(1,2,1)
        hold on

        y_values = [mean(all_glm(:,2)), mean(all_glm(:,3))];
        err_vals = [std(all_glm(:,2)), std(all_glm(:,3))];

        e = errorbar(1:2, y_values, err_vals, ...
            'k', 'LineWidth', 3, 'LineStyle', 'none');
        e.CapSize = 0;

        scatter(1:2, y_values, 180, 'k', 'filled')

        ylabel('Weights (a.u.)')
        xticks([])
        set(gca, 'FontSize', 15)

        %% --------------------------------------------------------
        %% Panel 2: Delta weight
        %% --------------------------------------------------------

        subplot(1,2,2)
        hold on

        delta_weight = abs(all_glm(:,2)) - abs(all_glm(:,3));

        e = errorbar(1, mean(delta_weight), std(delta_weight), ...
            'k', 'LineWidth', 3, 'LineStyle', 'none');
        e.CapSize = 0;

        scatter(1, mean(delta_weight), 180, 'k', 'filled')

        xlabel('\Delta weights (|win| - |lose|)')
        xticks([])
        set(gca, 'FontSize', 12)


    end
end
