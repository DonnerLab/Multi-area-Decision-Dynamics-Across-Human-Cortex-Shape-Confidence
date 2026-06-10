%% =======================================================================
% Figure S1D
% Correlation of contrast fluctuations across stimulus categories
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
    % STIMULUS CATEGORY
    %
    % Side = 1  --> Right > Left
    % Side = -1 --> Left > Right
    %% ===================================================================
    
    stimulusSide = double(mean(All_right,2) > mean(All_left,2));
    stimulusSide(stimulusSide == 0) = -1;
    
    
    %% ===================================================================
    % GENERATIVE MEAN FOR EACH TRIAL
    %% ===================================================================
    
    generativeMean = Real_data(:,4);
    
    
    %% ===================================================================
    % LEFT / RIGHT SAMPLE CONTRASTS
    %% ===================================================================
    
    rightSamples = All_right;
    leftSamples  = All_left;
    
    
    %% ===================================================================
    % REMOVE GENERATIVE MEAN
    %
    % Compute fluctuations around the generative mean
    %% ===================================================================
    
    rightSamples = ...
        rightSamples - (0.5 + generativeMean ./ 2);
    
    leftSamples = ...
        leftSamples - (0.5 - generativeMean ./ 2);
    
    
    %% ===================================================================
    % SPLIT TRIALS BY STIMULUS CATEGORY
    %% ===================================================================
    
    trials_RightGreater = find(stimulusSide == 1);
    trials_LeftGreater  = find(stimulusSide == -1);
    
    
    % Right > Left trials
    right_Rgreater = rightSamples(trials_RightGreater,:);
    left_Rgreater  = leftSamples(trials_RightGreater,:);
    
    
    % Left > Right trials
    right_Lgreater = rightSamples(trials_LeftGreater,:);
    left_Lgreater  = leftSamples(trials_LeftGreater,:);
    
    
    %% ===================================================================
    % CORRELATION OF CONTRAST FLUCTUATIONS
    %% ===================================================================
    
    % Category 1: Right > Left
    [corr_rightGreater, ~] = corr( ...
        left_Rgreater(:), ...
        right_Rgreater(:));
    
    
    % Category 2: Left > Right
    [corr_leftGreater, ~] = corr( ...
        left_Lgreater(:), ...
        right_Lgreater(:));
    
    
    %% ===================================================================
    % STORE RESULTS
    %% ===================================================================
    
    All_correlations(subject,:) = [ ...
        corr_rightGreater, ...
        corr_leftGreater];
    
end


%% =======================================================================
% FIGURE 
%% =======================================================================

figure
clf
hold on


%% =======================================================================
% GROUP MEAN BARS
%% =======================================================================

bar( ...
    1:2, ...
    mean(All_correlations), ...
    'facecolor', 'k', ...
    'facealpha', 0.5);


%% =======================================================================
% INDIVIDUAL SUBJECT POINTS
%% =======================================================================

for subject = 1:20
    
    for condition = 1:2
        
        scatter( ...
            condition + (rand(1)*2 - 1) * 0.05, ...
            All_correlations(subject,condition), ...
            100, ...
            'k', ...
            'filled', ...
            'markerfacealpha', 0.5);
        
    end
    
end


%% =======================================================================
% AXIS SETTINGS
%% =======================================================================

ylim([-1 1])

xticks([1 2])
xticklabels({'Left > Right','Right > Left'})
ylabel('Corr. Coeff.')
set(gca,'fontsize',15)

