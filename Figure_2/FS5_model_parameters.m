%% =======================================================================
% Model parameter summary across subjects
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
behaviorFolder = fullfile(projectFolder, 'data', 'behavior');
modelFolder    = fullfile(projectFolder, 'data', 'model', 'results');

%% =======================================================================
% Load model parameters
%% =======================================================================

for s = 1:20

    load(fullfile(modelFolder, ...
        ['Parameters_model_9_sbj_', num2str(s), '.mat']))

    Par_sbj = Parameters(s,:,2);

    % rescale noise parameter
    Par_sbj(2) = Par_sbj(2) ./ 20;

    All_pars(s,:) = Par_sbj;

    %% ---------------------------------------------------------------
    % Adaptation function
    %% ---------------------------------------------------------------

    tau_a = Par_sbj(3) * 15 + 0.0001;

    adaptation(s,:) = ...
        Par_sbj(4) + ...
        (1 - Par_sbj(4)) .* exp(-(0:9) ./ tau_a);

    adaptation(s,:) = adaptation(s,:) ./ sum(adaptation(s,:));

end

%% =======================================================================
% Plot average adaptation profile
%% =======================================================================

figure(100)

plot(mean(adaptation, 1), ...
    'color', 'k', ...
    'linewidth', 5)

hold on

xticks(1:10)
xlabel('Sample number')
ylabel('A.U.')
set(gca, 'fontsize', 15)
box off

%% =======================================================================
% Plot model parameters
%% =======================================================================

figure(1)

%% Lambda

subplot(2,7,2)

bar(1, mean(All_pars(:,1)), 'facecolor', 'r')
hold on

for s = 1:20
    scatter( ...
        1 + (rand(1)*2 - 1)*0.05, ...
        All_pars(s,1), ...
        40, ...
        'k', ...
        'filled', ...
        'markerfacealpha', 0.5)
end

ylabel('Lambda')
set(gca, 'fontsize', 15)

%% Adaptation decay

subplot(2,7,3)

bar(1, mean(All_pars(:,3)) * 15, 'facecolor', 'r')
hold on

for s = 1:20
    scatter( ...
        1 + (rand(1)*2 - 1)*0.05, ...
        All_pars(s,3) * 15, ...
        40, ...
        'k', ...
        'filled', ...
        'markerfacealpha', 0.5)
end

ylabel('Adaptation decay')
set(gca, 'fontsize', 15)

%% Decision bound offset

subplot(2,7,5)

bar(1, mean(All_pars(:,8)) * 15, 'facecolor', 'r')
hold on

for s = 1:20
    scatter( ...
        1 + (rand(1)*2 - 1)*0.05, ...
        All_pars(s,8) * 15, ...
        40, ...
        'k', ...
        'filled', ...
        'markerfacealpha', 0.5)
end

ylabel('Decision bound offset')
set(gca, 'fontsize', 15)

%% Bias

subplot(2,7,1)

bar(1, mean(All_pars(:,11)), 'facecolor', 'r')
hold on

for s = 1:20
    scatter( ...
        1 + (rand(1)*2 - 1)*0.05, ...
        All_pars(s,11), ...
        40, ...
        'k', ...
        'filled', ...
        'markerfacealpha', 0.5)
end

ylabel('Bias')
set(gca, 'fontsize', 15)

%% Noise

subplot(2,7,8)

bar(1, mean(All_pars(:,2)), 'facecolor', 'r')
hold on

for s = 1:20
    scatter( ...
        1 + (rand(1)*2 - 1)*0.05, ...
        All_pars(s,2), ...
        40, ...
        'k', ...
        'filled', ...
        'markerfacealpha', 0.5)
end

ylabel('Noise')
set(gca, 'fontsize', 15)

%% FFI weight

subplot(2,7,6)

bar(1, mean(All_pars(:,5)), 'facecolor', 'r')
hold on

for s = 1:20
    scatter( ...
        1 + (rand(1)*2 - 1)*0.05, ...
        All_pars(s,5), ...
        40, ...
        'k', ...
        'filled', ...
        'markerfacealpha', 0.5)
end

ylabel('FFI weight')
set(gca, 'fontsize', 15)

%% Accumulator balance weight

subplot(2,7,13)

bar(1, mean(All_pars(:,10)), 'facecolor', 'r')
hold on

for s = 1:20
    scatter( ...
        1 + (rand(1)*2 - 1)*0.05, ...
        All_pars(s,10), ...
        40, ...
        'k', ...
        'filled', ...
        'markerfacealpha', 0.5)
end

ylabel('Accumulator balance weight')
set(gca, 'fontsize', 15)

%% Confidence slope

subplot(2,7,7)

bar(1, mean(All_pars(:,7) * 10), 'facecolor', 'r')
hold on

for s = 1:20
    scatter( ...
        1 + (rand(1)*2 - 1)*0.05, ...
        All_pars(s,7) * 10, ...
        40, ...
        'k', ...
        'filled', ...
        'markerfacealpha', 0.5)
end

ylabel('Slope')
set(gca, 'fontsize', 15)

%% Confidence threshold

subplot(2,7,14)

bar(1, mean(All_pars(:,9) * 4), 'facecolor', 'r')
hold on

for s = 1:20
    scatter( ...
        1 + (rand(1)*2 - 1)*0.05, ...
        All_pars(s,9) * 4, ...
        40, ...
        'k', ...
        'filled', ...
        'markerfacealpha', 0.5)
end

ylabel('Threshold')
set(gca, 'fontsize', 15)

