%% =======================================================================
% Figure 1F-G: Plot PK analyses with confidence
%% =======================================================================

clear all
close all

% Folder where this script is located
scriptFolder = fileparts(mfilename('fullpath'));

% Go up two folders 
projectFolder = fileparts(fileparts(scriptFolder));

% Path to model results
modelFolder = fullfile(projectFolder, 'data', 'model', 'results');

% Path to functions
functionsFolder = fullfile(projectFolder, 'functions');
addpath(functionsFolder);


%% =======================================================================
% MAIN LOOP ACROSS SUBJECTS
%% =======================================================================

for subject = 1:20
    
    %% ===================================================================
    % LOAD MODEL + BEHAVIORAL DATA
    %% ===================================================================
    
    load(fullfile( ...
        modelFolder, ...
        ['All_data_model_9_sbj_' num2str(subject) '.mat']));
    
    
    %% ===================================================================
    % GENERATIVE MEAN FOR EACH TRIAL
    %% ===================================================================
    
    generativeMean_trial = Real_data(:,4);
    
    
    %% ===================================================================
    % MODEL INPUT
    %% ===================================================================
    
    simulationInput.t     = 10;
    simulationInput.right = All_right';
    simulationInput.left  = All_left';
    simulationInput.rt    = Real_data(:,3);
    
    
    %% ===================================================================
    % SIMULATED CONFIDENCE REPORTS
    %
    % Build multinomial samples from model confidence probabilities
    %% ===================================================================
    
    simulatedConfidenceChoices = [];
    
    for trial = 1:size(All_ch_model,1)
        
        simulatedConfidenceChoices(:,trial,:) = ...
            mnrnd(1, All_Conf(trial,:), 150);
        
    end
    
    simulatedConfidenceChoices = ...
        permute(simulatedConfidenceChoices, [3,2,1]);
    
    
    %% ===================================================================
    % LOOP OVER CHOICE AND SIDE
    %
    % choice:
    % 1 = right choice
    % 2 = left choice
    %
    % side:
    % 1 = left stimulus
    % 2 = right stimulus
    %% ===================================================================
    
    for choice = 1:2
        
        for side = 1:2
            
            %% ============================================================
            % SELECT SIDE-SPECIFIC CONTRASTS
            %% ============================================================
            
            if side == 1
                
                stimulusFluctuations = All_left;
                stimulusFluctuations_model = All_left;
                
                % Generative mean for left side
                generativeMean = 0.5 - generativeMean_trial ./ 2;
                
            else
                
                stimulusFluctuations = All_right;
                stimulusFluctuations_model = All_right;
                
                % Generative mean for right side
                generativeMean = 0.5 + generativeMean_trial ./ 2;
                
            end
            
            
            %% ============================================================
            % EXPAND VARIABLES FOR MODEL SIMULATIONS
            %% ============================================================
            
            stimulusFluctuations_model = ...
                repmat(stimulusFluctuations_model, [1,1,150]);
            stimulusFluctuations_model = ...
                permute(stimulusFluctuations_model, [2,1,3]);
            
            generativeMean_model = ...
                repmat(generativeMean, [1,1,150]);
            generativeMean_model = ...
                permute(generativeMean_model, [2,1,3]);
            
            stimulusFluctuations = ...
                repmat(stimulusFluctuations, [1,1,150]);
            stimulusFluctuations = ...
                permute(stimulusFluctuations, [2,1,3]);
            
            
            %% ============================================================
            % REAL CHOICE + CONFIDENCE
            %% ============================================================
            
            realChoice = Real_data(:,1);
            realConfidence = Real_data(:,2) - 1;
            
            realChoice = repmat(realChoice, [1,1,150]);
            realChoice = permute(realChoice, [2,1,3]);
            
            realConfidence = repmat(realConfidence, [1,1,150]);
            realConfidence = permute(realConfidence, [2,1,3]);
            
            
            %% ============================================================
            % CONFIDENCE MAPPING
            %% ============================================================
            
            confidenceMap = [2 1; 4 3];
            
            
            %% ============================================================
            % LOOP OVER LOW / HIGH CONFIDENCE
            %% ============================================================
            
            for confidenceLevel = 1:2
                
                if choice == 1
                    
                    % Real data
                    selectedTrials = ...
                        (realChoice == 1 & ...
                         realConfidence == confidenceLevel - 1);
                    
                    % Model
                    selectedTrials_model = ...
                        (simulatedConfidenceChoices( ...
                            confidenceMap(1, confidenceLevel), :, :) == 1);
                    
                else
                    
                    % Real data
                    selectedTrials = ...
                        (realChoice == 0 & ...
                         realConfidence == confidenceLevel - 1);
                    
                    % Model
                    selectedTrials_model = ...
                        (simulatedConfidenceChoices( ...
                            confidenceMap(2, confidenceLevel), :, :) == 1);
                    
                end
                
                
                %% ========================================================
                % COMPUTE CONTRAST FLUCTUATIONS
                %% ========================================================
                
                contrastFluctuation_model = mean( ...
                    (stimulusFluctuations_model(:,selectedTrials_model))' ...
                    - generativeMean_model(selectedTrials_model), ...
                    1);
                
                contrastFluctuation_real = mean( ...
                    (stimulusFluctuations(:,selectedTrials))' ...
                    - generativeMean_model(selectedTrials), ...
                    1);
                
                
                %% ========================================================
                % STORE PKS
                %% ========================================================
                
                PK_all(subject,:,1,choice,side,confidenceLevel) = ...
                    contrastFluctuation_real(1:10);
                
                PK_all_model(subject,:,1,choice,side,confidenceLevel) = ...
                    contrastFluctuation_model;
                
                PK_both(subject,:,1,confidenceLevel,side,choice) = ...
                    contrastFluctuation_real;
                
                PK_both_model(subject,:,1,confidenceLevel,side,choice) = ...
                    contrastFluctuation_model;
                
            end
        end
    end
end


%% =======================================================================
% COLORS
%% =======================================================================

lightBlue = [0 162 255] ./ 255;


%% =======================================================================
% EXTRACT PK FOR SELECTED / NON-SELECTED SIDE
%% =======================================================================

% REAL DATA
PK_selected = cat(3, ...
    PK_all(:,:,:,1,2,:), ...
    PK_all(:,:,:,2,1,:));

PK_selected = ...
    PK_selected(:,:,:,:,:,2) - PK_selected(:,:,:,:,:,1);

PK_selected = mean(PK_selected,3);


PK_nonSelected = cat(3, ...
    PK_all(:,:,:,1,1,:), ...
    PK_all(:,:,:,2,2,:));

PK_nonSelected = ...
    PK_nonSelected(:,:,:,:,:,2) - PK_nonSelected(:,:,:,:,:,1);

PK_nonSelected = mean(PK_nonSelected,3);


% MODEL
PK_selected_model = cat(3, ...
    PK_all_model(:,:,:,1,2,:), ...
    PK_all_model(:,:,:,2,1,:));

PK_selected_model = ...
    PK_selected_model(:,:,:,:,:,2) - PK_selected_model(:,:,:,:,:,1);

PK_selected_model = mean(PK_selected_model,3);


PK_nonSelected_model = cat(3, ...
    PK_all_model(:,:,:,1,1,:), ...
    PK_all_model(:,:,:,2,2,:));

PK_nonSelected_model = ...
    PK_nonSelected_model(:,:,:,:,:,2) - PK_nonSelected_model(:,:,:,:,:,1);

PK_nonSelected_model = mean(PK_nonSelected_model,3);


%% =======================================================================
% FIGURE 1 - DIFFERENCE BETWEEN SELECTED AND NON-SELECTED SIDE
%% =======================================================================

close all
figure(1)
hold on

AUC = [];
AUC_model = [];

for sample = 1:10
    
    AUC(:,sample) = ...
        abs(PK_selected(:,sample)) - abs(PK_nonSelected(:,sample));
    
    AUC_model(:,sample) = ...
        abs(PK_selected_model(:,sample)) - ...
        abs(PK_nonSelected_model(:,sample));
    
    
    e = errorbar( ...
        sample, ...
        nanmean(AUC(:,sample),1), ...
        nanstd(AUC(:,sample))./sqrt(20), ...
        'color','k', ...
        'linewidth',5, ...
        'linestyle','none');
    
    e.CapSize = 0;
    
    scatter( ...
        sample, ...
        nanmean(AUC(:,sample),1), ...
        80, ...
        'k', ...
        'filled', ...
        'MarkerEdgeColor',[1 1 1]);
    
    
    pValues(sample) = permutationTest( ...
        AUC(:,sample), ...
        zeros(20,1), ...
        100000);
    
end


significantP = double(fdr_bh(pValues,0.05));
significantP(significantP == 0) = nan;
significantP(significantP == 1) = 0.05;


shadedErrorBar( ...
    1:10, ...
    nanmean(AUC_model,1), ...
    nanstd(AUC_model,0,1)./sqrt(size(AUC_model,1)), ...
    'lineProps',{'.-','color','k'});

scatter( ...
    1:10, ...
    significantP ./ 3, ...
    80, ...
    'k', ...
    'square', ...
    'filled');

ylabel({'Excess contrast'; '(high confidence - low confidence)'})
xlabel('Sample number')

xticks(1:10)
xlim([0 11])

set(gca,'fontsize',15)


%% =======================================================================
% FIGURE 2 - PLOT EACH PSYCHOPHYSICAL KERNEL (PK)
%
% side = 1 --> selected side
% side = 2 --> non-selected side
%
% Low vs High confidence
% Real data + model predictions
%% =======================================================================

figure(2)
clf

for side = 1:2
    
    hold on
    
    
    %% ===================================================================
    % SIDE 1
    %% ===================================================================
    
    if side == 1
        
        % ---------------------------------------------------------------
        % REAL DATA
        % ---------------------------------------------------------------
        
        PK_selectedSide = cat( ...
            5, ...
            PK_both(:,:,1,:,1,1), ...
            PK_both(:,:,1,:,2,2));
        
        
        % Low confidence
        e = errorbar( ...
            1:10, ...
            nanmean(PK_selectedSide(:,:,1,1),1), ...
            nanstd(PK_selectedSide(:,:,1,1))./sqrt(20), ...
            'color', [0 0.4470 0.7410], ...
            'linewidth', 5, ...
            'linestyle', 'none');
        
        e.CapSize = 0;
        
        scatter( ...
            1:10, ...
            nanmean(PK_selectedSide(:,:,1,1),1), ...
            80, ...
            [0 0.4470 0.7410], ...
            'filled', ...
            'MarkerEdgeColor', [1 1 1]);
        
        
        % High confidence
        e = errorbar( ...
            1:10, ...
            nanmean(PK_selectedSide(:,:,1,2),1), ...
            nanstd(PK_selectedSide(:,:,1,2))./sqrt(20), ...
            'color', lightBlue, ...
            'linewidth', 5, ...
            'linestyle', 'none');
        
        e.CapSize = 0;
        
        scatter( ...
            1:10, ...
            nanmean(PK_selectedSide(:,:,1,2),1), ...
            80, ...
            lightBlue, ...
            'filled', ...
            'MarkerEdgeColor', [1 1 1], ...
            'markerfacealpha', 0.5);
        
        
        % ---------------------------------------------------------------
        % MODEL
        % ---------------------------------------------------------------
        
        PK_selectedSide_model = cat( ...
            5, ...
            PK_both_model(:,:,1,:,1,1), ...
            PK_both_model(:,:,1,:,2,2));
        
        
        % Low confidence
        shadedErrorBar( ...
            1:10, ...
            nanmean( ...
                nanmean(PK_selectedSide_model(:,:,1,1,:),5),1), ...
            nanstd( ...
                nanmean(PK_selectedSide_model(:,:,1,1,:),5),0,1) ...
                ./ sqrt(size(PK_selectedSide_model,1)), ...
            'lineProps', {'-','color',[0 0.4470 0.7410]});
        
        
        % High confidence
        shadedErrorBar( ...
            1:10, ...
            nanmean( ...
                nanmean(PK_selectedSide_model(:,:,1,2,:),5),1), ...
            nanstd( ...
                nanmean(PK_selectedSide_model(:,:,1,2,:),5),0,1) ...
                ./ sqrt(size(PK_selectedSide_model,1)), ...
            'lineProps', {'-','color', lightBlue});
        
        
    %% ===================================================================
    % SIDE 2
    %% ===================================================================
    
    else
        
        % ---------------------------------------------------------------
        % REAL DATA
        % ---------------------------------------------------------------
        
        PK_nonSelectedSide = cat( ...
            5, ...
            PK_both(:,:,1,:,1,2), ...
            PK_both(:,:,1,:,2,1));
        
        
        % Low confidence
        e = errorbar( ...
            1:10, ...
            nanmean(PK_nonSelectedSide(:,:,1,1),1), ...
            nanstd(PK_nonSelectedSide(:,:,1,1))./sqrt(20), ...
            'color', 'k', ...
            'linewidth', 5, ...
            'linestyle', 'none');
        
        e.CapSize = 0;
        
        scatter( ...
            1:10, ...
            nanmean(PK_nonSelectedSide(:,:,1,1),1), ...
            80, ...
            'k', ...
            'filled', ...
            'MarkerEdgeColor', [1 1 1]);
        
        
        % High confidence
        e = errorbar( ...
            1:10, ...
            nanmean(PK_nonSelectedSide(:,:,1,2),1), ...
            nanstd(PK_nonSelectedSide(:,:,1,2))./sqrt(20), ...
            'color', [0.5 0.5 0.5], ...
            'linewidth', 5, ...
            'linestyle', 'none');
        
        e.CapSize = 0;
        
        scatter( ...
            1:10, ...
            nanmean(PK_nonSelectedSide(:,:,1,2),1), ...
            80, ...
            [0.5 0.5 0.5], ...
            'filled', ...
            'MarkerEdgeColor', [1 1 1]);
        
        
        % ---------------------------------------------------------------
        % MODEL
        % ---------------------------------------------------------------
        
        PK_nonSelectedSide_model = cat( ...
            5, ...
            PK_both_model(:,:,1,:,1,2), ...
            PK_both_model(:,:,1,:,2,1));
        
        
        % Low confidence
        shadedErrorBar( ...
            1:10, ...
            nanmean( ...
                nanmean(PK_nonSelectedSide_model(:,:,1,1,:),5),1), ...
            nanstd( ...
                nanmean(PK_nonSelectedSide_model(:,:,1,1,:),5),0,1) ...
                ./ sqrt(size(PK_nonSelectedSide_model,1)), ...
            'lineProps', {'-','color','k'});
        
        
        % High confidence
        shadedErrorBar( ...
            1:10, ...
            nanmean( ...
                nanmean(PK_nonSelectedSide_model(:,:,1,2,:),5),1), ...
            nanstd( ...
                nanmean(PK_nonSelectedSide_model(:,:,1,2,:),5),0,1) ...
                ./ sqrt(size(PK_nonSelectedSide_model,1)), ...
            'lineProps', {'-','color',[0.5 0.5 0.5]});
        
    end
    
    
    %% ===================================================================
    % AXIS SETTINGS
    %% ===================================================================
    
    ylabel('Excess contrast')
    xlabel('Sample N')
    
    xticks(1:10)
    xlim([0 11])
    
    set(gca,'fontsize',15)
    
end