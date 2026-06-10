%% Example neural DV trajectories
clear;
close all;
clc;


%% =======================================================================
% Paths
%% =======================================================================
% Folder where this script is located
scriptFolder = fileparts(mfilename('fullpath'));

% Go up two folders 
projectFolder = fileparts(fileparts(scriptFolder));

% Main folders
data_path      = fullfile(projectFolder, 'data', 'decoder_output', 'hemi');
decoderFolder  = fullfile(projectFolder, 'data', 'brain_map_decoder');
functionsFolder = fullfile(projectFolder, 'functions');
giftiFolder     = fullfile(functionsFolder, 'gifti-1.6');

% Add required paths
addpath(data_path)
addpath(decoderFolder)
addpath(functionsFolder)
addpath(giftiFolder)

% Move to data folder
cd(data_path)




cd(data_path);

%% Subject/session info
SubjectID_2 = { ...
    'S01','S04','S06','S07','S12','S15','S17','S18', ...
    'S11','S05','S19','S20','S22','S10','S02','S08', ...
    'S24','S21','S23','S25'};

SubjectID = [1,4,6,7,12,15,17,18,11,5,19,20,22,10,2,8,24,21,23,25];

sessions = [ ...
    11,2,3,4,5,6;
    11,2,3,4,5,6;
    1,2,3,4,5,6;
    11,2,3,4,5,6;
    11,2,3,4,5,6;
    11,2,3,4,5,6;
    11,2,3,4,5,6;
    11,2,3,4,5,6;
    11,2,3,4,5,6;
    11,2,3,4,5,6;
    1,2,3,4,5,6;
    1,2,3,4,5,6;
    1,2,3,4,5,6;
    11,2,3,4,5,6;
    11,2,3,4,5,6;
    1,2,3,4,5,6;
    1,2,3,4,5,6;
    1,2,3,4,5,6;
    1,2,3,4,5,6;
    1,2,3,4,5,6];

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

%% Hemisphere labels
hemis = {'l','r'};

%% Example subject/session only
% Original code uses:
%   side_n = 1
%   s = 1
%   ses = 1

subject_idx = 1;
session_idx = 1;

%% Build file names
file_left = [ ...
    SubjectID_2{subject_idx}, ...
    '_choice_', hemis{1}, ...
    '_session', num2str(session_idx), ...
    '_PMd_M1_decoder_output_100hz.mat'];

file_right = [ ...
    SubjectID_2{subject_idx}, ...
    '_choice_', hemis{2}, ...
    '_session', num2str(session_idx), ...
    '_PMd_M1_decoder_output_100hz.mat'];

if ~isfile(file_left) || ~isfile(file_right)
    error('Decoder files not found for selected subject/session.');
end

%% Load decoder output: left hemisphere
load(file_left);

Perf_L = Prob_Performance;
Perf_L = permute(Perf_L, [3, 2, 1]);

%% Load decoder output: right hemisphere
load(file_right);

Perf_R = Prob_Performance;
Perf_R = permute(Perf_R, [3, 2, 1]);

%% Trial-wise variables
Evidence = abs( ...
    mean(All_samples(:,1:10), 2) - ...
    mean(All_samples(:,11:20), 2));

Q = quantile(Evidence(:,1), 2);

Side = -sign( ...
    mean(All_samples(:,1:10), 2) - ...
    mean(All_samples(:,11:20), 2));

Correctness = (All_Choices == Side);

[N, edges, id_evidence] = histcounts( ...
    Evidence(:,1), ...
    [-inf -flip(Q) 0 Q inf]);

%% Trial selection examples

% Congruent:
% choice left + high confidence
% right hemi supports left choice
% left hemi supports right choice

id_congr_Rh_Cl = find( ...
    mean(Perf_R(:,141:151),2) > 0 & ...
    mean(Perf_L(:,141:151),2) < 0 & ...
    All_Choices == -1 & ...
    Evidence > 0.07 & ...
    All_Conf == 2);

% Incongruent:
% choice left + low confidence
% both hemis support same side

id_inc_Rh_Cl = find( ...
    mean(Perf_R(:,141:151),2) > 0 & ...
    mean(Perf_L(:,141:151),2) > 0 & ...
    All_Choices == -1 & ...
    Evidence > 0.07 & ...
    All_Conf == 1);

%% Manually selected example trials
id_congr_CL   = 234;
id_incongr_CL = 462;

Rt_congr_CL   = All_RT(id_congr_CL);
Rt_incongr_CL = All_RT(id_incongr_CL);

%% Smooth DV during stimulus presentation
stim_window = 51:151;
smooth_k = 25;

R_h_congr_CL = movmean( ...
    mean(Perf_R(id_congr_CL, stim_window), 1), ...
    smooth_k);

L_h_congr_CL = movmean( ...
    mean(Perf_L(id_congr_CL, stim_window), 1), ...
    smooth_k);

R_h_incongr_CL = movmean( ...
    mean(Perf_R(id_incongr_CL, stim_window), 1), ...
    smooth_k);

L_h_incongr_CL = movmean( ...
    mean(Perf_L(id_incongr_CL, stim_window), 1), ...
    smooth_k);

%% Plot trajectory in neural state space
figure(99);
clf;

subplot(1,2,1);

% Congruent trial
plot( ...
    R_h_congr_CL, ...
    L_h_congr_CL, ...
    'color', 'k', ...
    'linewidth', 2);
hold on;

scatter( ...
    R_h_congr_CL(end), ...
    L_h_congr_CL(end), ...
    50, ...
    'k', ...
    'filled');

% Incongruent trial
plot( ...
    R_h_incongr_CL, ...
    L_h_incongr_CL, ...
    'color', [0.5 0.5 0.5], ...
    'linewidth', 2);

scatter( ...
    R_h_incongr_CL(end), ...
    L_h_incongr_CL(end), ...
    50, ...
    [0.5 0.5 0.5], ...
    'filled');

%% Formatting
xlabel('Acc. evidence for left choice');
ylabel('Acc. evidence for right choice');

xline(0);
yline(0);

axis square;
xlim([-0.5 0.5]);
ylim([-0.5 0.5]);

set(gca, 'fontsize', 15);
box off;