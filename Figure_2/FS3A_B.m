%% Clean Decoder Performance Analysis
% Stimulus-aligned + response-aligned decoder performance
% Single-trial slope analysis and evidence-strength correlations

clear; close all; clc;

%% Paths and subject IDs

% path generation
scriptFolder  = fileparts(mfilename('fullpath'));
projectFolder = fileparts(fileparts(scriptFolder));

base_path = fullfile(projectFolder, 'data', 'decoder_output', 'bilateral');
response_path = fullfile(base_path, 'response');
stimulus_path = fullfile(base_path, 'stimulus');

% Required helper functions
functionsPath = fullfile(projectFolder, 'functions');
addpath(functionsPath);


subject_ids = { ...
    'S01','S04','S06','S07','S12','S15','S17','S18', ...
    'S11','S05','S19','S20','S22','S10','S02','S08', ...
    'S24','S21','S23','S25'};

n_subjects = 20;
n_sessions = 6;

%% ================================================================
%% PART 1: RESPONSE-ALIGNED DECODER PERFORMANCE
%% ================================================================

response_path = fullfile(base_path, 'response');
cd(response_path)

response_output = nan(n_subjects, n_sessions, 401);

for s = 1:n_subjects
    for ses = 1:n_sessions

        file_name = sprintf('%s_choice_both_session%d_PCA_PMd_M1_resp.mat', ...
            subject_ids{s}, ses);

        if isfile(file_name)
            load(file_name, 'Decoder_Performance');
            response_output(s, ses, :) = Decoder_Performance;
        end
    end
end

response_output(response_output == 0) = NaN;

%% Plot response-aligned decoder performance
figure;
subplot(1,5,5)

response_time = -0.5:0.0025:0.5;
response_mean = reshape(nanmean(nanmean(response_output,1),2), 1, 401);
response_sem  = reshape( ...
    nanstd(nanmean(response_output,2),0,1), 1, 401) ./ sqrt(n_subjects);

shadedErrorBar(response_time, response_mean, response_sem, ...
    'lineProps', {'-','color','k'});

xlim([-0.3 0.2])
ylim([0.48 0.8])
xline(0, '--', 'LineWidth', 3, 'Color', [0.5 0.5 0.5])

xlabel('Time from response (s)')
set(gca, 'FontSize', 15)


%% ================================================================
%% PART 2: STIMULUS-ALIGNED DECODER PERFORMANCE
%% ================================================================

stimulus_path = fullfile(base_path, 'stimulus');
cd(stimulus_path)

stim_output = nan(n_subjects, n_sessions, 201);

for s = 1:n_subjects
    for ses = 1:n_sessions

        file_name = sprintf('%s_choice_both_session%d_PCA_PMd_M1_100Hz.mat', ...
            subject_ids{s}, ses);

        if isfile(file_name)
            load(file_name, 'Decoder_Performance');
            stim_output(s, ses, :) = Decoder_Performance;
        end
    end
end

stim_output(stim_output == 0) = NaN;

%% Permutation test against chance (0.5)
all_subject_mean = reshape(mean(stim_output,2), n_subjects, 201);
p_values = nan(1,201);

for t = 1:201
    p_values(t) = permutationTest( ...
        all_subject_mean(:,t) - 0.5, ...
        zeros(n_subjects,1), ...
        10000);
end

significant_points = fdr_bh(p_values, 0.05);

%% Plot stimulus-aligned decoder performance
figure(1)
subplot(1,5,1:4)

stim_time = -0.5:0.01:1.5;
stim_mean = reshape(nanmean(nanmean(stim_output,1),2), 1, 201);
stim_sem  = reshape( ...
    nanstd(nanmean(stim_output,2),0,1), 1, 201) ./ sqrt(n_subjects);

shadedErrorBar(stim_time, stim_mean, stim_sem, ...
    'lineProps', {'-','color','k'});
hold on

scatter(stim_time, significant_points .* 0.7, ...
    150, 'k', 'filled', 'square')

xline(0, '--', 'LineWidth', 3, 'Color', [0.5 0.5 0.5])
xline(1, '--', 'LineWidth', 3, 'Color', [0.5 0.5 0.5])

for i = 1:9
    xline(0.1 + 0.1*(i-1), '--', ...
        'LineWidth', 2, ...
        'Color', [0.5 0.5 0.5]);
end

xlim([-0.5 1.3])
ylim([0.48 0.8])

xlabel('Time from stimulus onset (s)')
set(gca, 'FontSize', 15)


%% ================================================================
%% PART 3: DECODER SLOPE VS EVIDENCE STRENGTH
%% ================================================================

slope_corr = nan(n_subjects, n_sessions);
end_corr   = nan(n_subjects, n_sessions);
quantile_output = nan(n_subjects, n_sessions, 201, 3);

for s = 1:n_subjects
    for ses = 1:n_sessions

        perf_file = sprintf('%s_choice_both_session%d_PCA_PMd_M1_100Hz.mat', ...
            subject_ids{s}, ses);

        decoder_file = sprintf('%s_choice_both_session%d_PMd_M1_decoder_output_100hz.mat', ...
            subject_ids{s}, ses);

        if ~(isfile(perf_file) && isfile(decoder_file))
            continue
        end

        load(perf_file)
        load(decoder_file)

        %% Evidence strength
        evidence = abs(mean(All_samples(:,1:10) - All_samples(:,11:20), 2));
        [~,~,evidence_bin] = histcounts(evidence, [-inf 0.06 0.12 inf]);

        %% Decoder probability
        neural_dv = 1 ./ (1 + exp(Prob_Performance));
        binary_dv = neural_dv > 0.5;

        choices = repmat(All_Choices, 1, 201);
        choices(choices == -1) = 0;

        predicted = permute(binary_dv, [3 2 1]);
        correct_decoder = (predicted == choices);

        %% Trial-wise decoder confidence
        perf = nan(size(Prob_Performance));

        for tr = 1:size(choices,1)
            if choices(tr) == 1
                perf(:,:,tr) = 1 ./ (1 + exp(Prob_Performance(:,:,tr)));
            else
                perf(:,:,tr) = 1 - 1 ./ (1 + exp(Prob_Performance(:,:,tr)));
            end
        end

        perf_trials = permute(perf, [3 2 1]);

        %% Final decoder performance (last 100 ms)
        decoder_end = mean(perf_trials(:,141:151), 2);

        end_corr(s,ses) = corr( ...
            zscore(decoder_end - 0.5), ...
            zscore(evidence), ...
            'type', 'spearman');

        %% Decoder slope (200 ms to 1 s)
        trial_slope = nan(size(perf_trials,1),1);

        for tr = 1:size(perf_trials,1)
            fit_par = polyfit( ...
                (0.2:0.01:1)', ...
                (perf_trials(tr,71:151) - 0.5)', ...
                1);

            trial_slope(tr) = fit_par(1);
        end

        slope_corr(s,ses) = corr( ...
            zscore(trial_slope), ...
            zscore(evidence), ...
            'type', 'spearman');

        %% Performance by evidence quantile
        for q = 1:3
            quantile_output(s,ses,:,q) = ...
                mean(correct_decoder(evidence_bin == q,:), 1);
        end
    end
end


%% ================================================================
%% PART 4: PLOT SLOPE-EVIDENCE CORRELATION
%% ================================================================

figure;
mean_slope_corr = mean(slope_corr, 2);

bar(1, mean(mean_slope_corr), ...
    'FaceColor', 'k', ...
    'FaceAlpha', 0.5)
hold on

for s = 1:n_subjects
    scatter(1, mean_slope_corr(s), ...
        100, 'k', 'filled', ...
        'MarkerFaceAlpha', 0.5)
end

xticks([])
ylabel('Correlation Coefficient')
set(gca, 'FontSize', 15)


%% ================================================================
%% PART 5: PLOT EVIDENCE QUANTILES
%% ================================================================

quantile_output(quantile_output == 0) = NaN;

figure(1)
subplot(1,5,1:4)
hold on

colors = flip(colormap(brewermap(6,'*RdBu')), 1);

for q = 1:3

    curve = reshape( ...
        nanmean(nanmean(quantile_output(:,:,:,q),1),2), ...
        1, 201);

    curve = smoothdata(curve, 'gaussian', [10]);

    plot(stim_time, curve, ...
        'Color', colors(q,:), ...
        'LineWidth', 3)

    fit_par = polyfit( ...
        (0.2:0.01:1)', ...
        reshape( ...
            nanmean(nanmean(quantile_output(:,:,71:151,q),1),2), ...
            1, 81)', ...
        1);

    fit_line = smoothdata( ...
        (0.2:0.01:1) .* fit_par(1) + fit_par(2), ...
        'gaussian', [10]);

    plot((0.2:0.01:1), fit_line, ...
        'Color', colors(q,:), ...
        'LineWidth', 3)
end

xline(0, '--', 'LineWidth', 3, 'Color', [0.5 0.5 0.5])
xline(1, '--', 'LineWidth', 3, 'Color', [0.5 0.5 0.5])

for i = 1:9
    xline(0.1 + 0.1*(i-1), '--', ...
        'LineWidth', 2, ...
        'Color', [0.5 0.5 0.5]);
end

xlim([-0.5 1.3])
xlabel('Time from stimulus onset (s)')
set(gca, 'FontSize', 15)
