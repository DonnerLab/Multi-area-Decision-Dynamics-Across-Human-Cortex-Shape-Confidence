%% =======================================================================
% Figure 2F
%
%%===========================================================

clear all
close all


%% =======================================================================
% PATHS
%% =======================================================================

% Folder where this script is located
scriptFolder = fileparts(mfilename('fullpath'));
if isempty(scriptFolder), scriptFolder = pwd; end

% Go up one folder (repository root)
projectFolder = fileparts(scriptFolder);

% Single-hemisphere decoder output (PMd/M1, 100 Hz)
decoderFolder = fullfile(projectFolder, 'data', 'decoder_output', 'hemi');


%% =======================================================================
% EXAMPLE SUBJECT AND SESSION
%% =======================================================================

subjectID = 'S01';
session   = 1;

% Left hemisphere decoder: log odds > 0 = evidence for a right choice
load(fullfile(decoderFolder, [subjectID '_choice_l_session' num2str(session) ...
    '_PMd_M1_decoder_output_100hz.mat']));
DV_right = permute(Prob_Performance, [3 2 1]);   % trials x time

% Right hemisphere decoder: log odds > 0 = evidence for a left choice
load(fullfile(decoderFolder, [subjectID '_choice_r_session' num2str(session) ...
    '_PMd_M1_decoder_output_100hz.mat']));
DV_left = permute(Prob_Performance, [3 2 1]);    % trials x time


%% =======================================================================
% CANDIDATE TRIALS (left choice, evidence strength > 0.07)
%
% DVs averaged in the last 100 ms of stimulus viewing (samples 141:151)
% All_samples: trials x 20 contrasts (1:10 left, 11:20 right)
% All_Choices: -1 = left; All_Conf: 1 = low, 2 = high
%% =======================================================================

evidenceStrength = abs(mean(All_samples(:,1:10),2) - mean(All_samples(:,11:20),2));

DV_left_end  = mean(DV_left(:,141:151),2);
DV_right_end = mean(DV_right(:,141:151),2);

% Consistent: DV_left > 0 and DV_right < 0, high confidence
consistentTrials = find(DV_left_end > 0 & DV_right_end < 0 & ...
    All_Choices == -1 & evidenceStrength > 0.07 & All_Conf == 2);

% Inconsistent: DV_left > 0 and DV_right > 0, low confidence
inconsistentTrials = find(DV_left_end > 0 & DV_right_end > 0 & ...
    All_Choices == -1 & evidenceStrength > 0.07 & All_Conf == 1);

% Example trials selected among the candidates
consistentTrial   = 234;
inconsistentTrial = 462;


%% =======================================================================
% SMOOTHED DV TRAJECTORIES DURING STIMULUS VIEWING (samples 51:151)
%% =======================================================================

stimulusWindow = 51:151;
smoothing      = 25;   % moving-average window (samples)

trajectory_left_consistent    = movmean(DV_left(consistentTrial,stimulusWindow),    smoothing);
trajectory_right_consistent   = movmean(DV_right(consistentTrial,stimulusWindow),   smoothing);
trajectory_left_inconsistent  = movmean(DV_left(inconsistentTrial,stimulusWindow),  smoothing);
trajectory_right_inconsistent = movmean(DV_right(inconsistentTrial,stimulusWindow), smoothing);


%% =======================================================================
% FIGURE (Figure 2F)
%
% RTs of the two trials: All_RT(consistentTrial), All_RT(inconsistentTrial)
%% =======================================================================

figure
hold on

plot(trajectory_left_consistent, trajectory_right_consistent, 'color', 'k', 'linewidth', 2);
scatter(trajectory_left_consistent(end), trajectory_right_consistent(end), 50, 'k', 'filled');

plot(trajectory_left_inconsistent, trajectory_right_inconsistent, 'color', [0.5 0.5 0.5], 'linewidth', 2);
scatter(trajectory_left_inconsistent(end), trajectory_right_inconsistent(end), 50, [0.5 0.5 0.5], 'filled');

title(sprintf('RT: %g (black), %g (gray)', ...
    All_RT(consistentTrial), All_RT(inconsistentTrial)))

xline(0);
yline(0);
axis square
xlim([-0.5 0.5])
ylim([-0.5 0.5])
xlabel('DV\_Left (log odds)')
ylabel('DV\_Right (log odds)')
box off
set(gca, 'fontsize', 15);
