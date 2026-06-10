%% =======================================================================
% Figure S1E-F
% Psychometric curve + Psychophysical Kernel (PK)
%% =======================================================================

clear all
close all


%% =======================================================================
% SUBJECT INFORMATION
%% =======================================================================

Subjects = [1,2,4,5,7,8,11,12,15,17,18,10,6,19,20,21,22,23,24,25];

Drugs = [ ...
    3 1 1 2 3 2;
    3 1 3 2 1 2;
    1 3 2 1 3 2;
    2 1 1 3 2 3;
    1 3 2 3 1 2;
    3 1 1 2 2 3;
    1 2 2 3 3 1;
    3 1 2 1 3 2;
    1 3 3 2 2 1;
    2 1 3 2 3 1;
    1 2 3 2 1 3;
    1 2 2 3 1 3;
    2 1 3 1 2 3;
    2 3 1 2 1 3;
    1 2 3 1 3 2;
    2 1 1 2 3 3;
    3 1 2 2 1 3;
    2 3 2 1 1 3;
    3 1 1 2 3 2;
    3 2 1 3 1 2];

Extra_sessions = [ ...
    1,2,4,5,7,0,11,12,15,17,18,10,0,0,0,0,0,0,0,0,0];


%% =======================================================================
% PATHS
%% =======================================================================

% Folder where this script is located
scriptFolder = fileparts(mfilename('fullpath'));

% Go up two folders 
projectFolder = fileparts(fileparts(scriptFolder));

% Path to model results
modelFolder = fullfile(projectFolder, 'data', 'model', 'results');


%% =======================================================================
% MAIN LOOP ACROSS SUBJECTS
%% =======================================================================

for subject = 1:20
    
    %% ===================================================================
    % LOAD SUBJECT DATA
    %% ===================================================================
    
    load(fullfile( ...
        modelFolder, ...
        ['All_data_model_9_sbj_' num2str(subject) '.mat']));
    
    
    %% ===================================================================
    % MODEL CHOICES
    %% ===================================================================
    
    modelChoices = All_ch_model(:,1);
    
    
    %% ===================================================================
    % REAL CHOICES + EVIDENCE
    %% ===================================================================
    
    % Mean contrast difference (Right - Left)
    evidence = mean(All_right,2) - mean(All_left,2);
    
    % Real participant choices
    realChoices = Real_data(:,1);
    
    
    %% ===================================================================
    % PSYCHOMETRIC FUNCTION
    %
    % Fit psychometric curve with lapse parameters
    %% ===================================================================
    
    stimulus = @(b,x) b(1) .* x(:,1);
    
    psychometricFunction = @(b,x) ...
        (b(14) + ...
        (1 - b(15) - b(14)) .* ...
        (1 ./ (1 + exp(-(stimulus(b,x))))));
    
    
    % Initial values
    x0 = [10 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    
    % Bounds
    lb = [-5 -2 -2 -2 -2 -2 -2 -2 -2 -2 -2 -2 -2 0 0];
    ub = [100 2 2 2 2 2 2 2 2 2 2 2 2 0.5 0.5];
    
    
    %% ===================================================================
    % REMOVE NaNs
    %% ===================================================================
    
    x = evidence;
    y = realChoices';
    
    nanTrials = unique([ ...
        find(sum(isnan(x),2) > 0); ...
        find(isnan(y) > 0)']);
    
    x(nanTrials,:) = [];
    y(nanTrials) = [];
    modelChoices(nanTrials,:) = [];
    
    
    %% ===================================================================
    % FIT PSYCHOMETRIC CURVE
    %% ===================================================================
    
    fittedParameters = lsqcurvefit( ...
        psychometricFunction, ...
        x0, ...
        x, ...
        y', ...
        lb, ...
        ub);
    
    fittedChoices = psychometricFunction(fittedParameters, x);
    
    
    %% ===================================================================
    % BIN BY EVIDENCE QUANTILES
    %% ===================================================================
    
    quantiles = round(quantile(abs(x(:,1)),4),2);
    
    [~, edges, bins] = histcounts( ...
        x(:,1), ...
        [flip(-quantiles) quantiles]);
    
    
    for bin = 1:length(edges)-1
        
        % Real choices
        P_choice(subject,bin,1) = ...
            nanmean(y(bins == bin));
        
        % Psychometric fit
        P_choice_fit(subject,bin,1) = ...
            nanmean(fittedChoices(bins == bin));
        
        % Model choices
        P_choice_model(subject,bin,1) = ...
            nanmean(modelChoices(bins == bin));
        
        % Bin center
        binCenters(bin,:) = mean(edges(bin:bin+1));
        
    end
    
    
    %% ===================================================================
    % PSYCHOPHYSICAL KERNEL (PK)
    %% ===================================================================
    
    contrastDifference = All_right - All_left;
    
    stimulusSide = sign(x);
    
    
    % ---------------------------------------------------------------
    % Trials where Right > Left
    % ---------------------------------------------------------------
    
    for sample = 1:10
        
        [~,~,~,AUC] = perfcurve( ...
            Real_data(stimulusSide > 0,1), ...
            contrastDifference(stimulusSide > 0,sample), ...
            1);
        
        AUC_right(:,sample) = AUC;
        
    end
    
    
    % ---------------------------------------------------------------
    % Trials where Left > Right
    % ---------------------------------------------------------------
    
    for sample = 1:10
        
        [~,~,~,AUC] = perfcurve( ...
            Real_data(stimulusSide < 0,1), ...
            contrastDifference(stimulusSide < 0,sample), ...
            1);
        
        AUC_left(:,sample) = AUC;
        
    end
    
    
    % Final PK measure
    AUC_all(subject,:) = mean([AUC_left; AUC_right],1) - 0.5;
    
end


%% =======================================================================
% FIGURE 1 PSYCHOMETRIC CURVE
%% =======================================================================

figure(1)
clf
hold on


% ---------------------------------------------------------------
% Real data
% ---------------------------------------------------------------

for bin = 1:size(P_choice_fit,2)
    
    errorbar( ...
        binCenters(bin), ...
        nanmean(P_choice(:,bin,1),1), ...
        nanstd(P_choice(:,bin,1),1), ...
        'linewidth', 3, ...
        'linestyle', 'none', ...
        'color', 'k');
    
    scatter( ...
        binCenters(bin), ...
        nanmean(P_choice(:,bin,1),1), ...
        100, ...
        'markeredgecolor', 'k', ...
        'markerfacecolor', 'k');
    
end


% ---------------------------------------------------------------
% Model choices
% ---------------------------------------------------------------

for subject = 1:20
    
    plot( ...
        binCenters, ...
        nanmean(P_choice_model(subject,:,1),1), ...
        'color', [0.5 0.5 0.5], ...
        'linewidth', 2);
    
end


% Group mean
plot( ...
    binCenters, ...
    nanmean(P_choice_model(:,:,1),1), ...
    'linewidth', 5, ...
    'color', 'k');


box off

ylabel('Probability right choice')
xlabel('\Delta Contrast (R-L)')

xticks(binCenters)

set(gca,'fontsize',15)


%% =======================================================================
% FIGURE 2  PSYCHOPHYSICAL KERNEL (PK)
%% =======================================================================

figure(2)
clf
hold on


for subject = 1:20
    
    plot( ...
        1:10, ...
        nanmean(AUC_all(subject,:,1),1), ...
        'color', [0.5 0.5 0.5], ...
        'linewidth', 2);
    
end


% Group mean
plot( ...
    1:10, ...
    nanmean(AUC_all(:,:,1),1), ...
    'linewidth', 5, ...
    'color', 'k');


box off

ylabel('ROC index - 0.5')
xlabel('Sample number')

set(gca,'fontsize',15)