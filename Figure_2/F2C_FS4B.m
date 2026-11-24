%% =======================================================================
% Figure 2C (PMd), Figure S4B (M1)
%
% Temporal cross-decoding (generalization) matrix of the bilateral choice
% decoder
%% =======================================================================

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

% Generalization matrices, one file per subject and area group
matrixFolder = fullfile(projectFolder, 'data', 'generalization_matrix');


%% =======================================================================
% SETTINGS
%% =======================================================================

groups      = {'PMd', 'M1'};
groupTitles = {'PMd (6a + 6d)', 'M1 (area 4)'};

% Color range of the decoder accuracy
accuracyLimits = [0.42 0.58];

% Time axis of the decoded samples: 400 Hz, downsampled by 4 (10 ms
% bins), 201 samples from -0.5 to 1.5 s from stimulus onset
time = -0.5 + (0:200) * 4/400;

% Diverging red-blue colormap (blue < 0.5 < red)
redBlue = [  5  48  97;  33 102 172;  67 147 195; 146 197 222; 209 229 240; ...
           247 247 247; 253 219 199; 244 165 130; 214  96  77; 178  24  43; 103   0  31] / 255;
colorMap = interp1(linspace(0,1,size(redBlue,1)), redBlue, linspace(0,1,256));


%% =======================================================================
% FIGURE (Figure 2C: PMd, Figure S4B: M1)
%% =======================================================================

figure('Color','w', 'Position', [100 100 1050 430]);

for group = 1:numel(groups)

    % Load the matrix of every file in the folder (training time x testing time)
    files = dir(fullfile(matrixFolder, ['S*_' groups{group} '_CHOICEdecoder_GENERALIZATION.mat']));

    allMatrices = [];

    for file = 1:numel(files)
        load(fullfile(matrixFolder, files(file).name), 'accuracy_gen');
        allMatrices = cat(3, allMatrices, accuracy_gen);
    end

    % Group average
    groupMatrix = mean(allMatrices, 3, 'omitnan');

    % Rows = training time, columns = testing time
    ax = subplot(1, numel(groups), group);
    imagesc(time, time, groupMatrix);
    axis xy square
    colormap(ax, colorMap);
    caxis(accuracyLimits);
    hold on

    % Diagonal and stimulus onset
    plot(time([1 end]), time([1 end]), 'k-', 'LineWidth', 0.5);
    xline(0, 'k--', 'LineWidth', 1);
    yline(0, 'k--', 'LineWidth', 1);

    xlabel('Testing time from stimulus onset (s)');
    ylabel('Training time from stimulus onset (s)');
    title(sprintf('%s   (n = %d)', groupTitles{group}, size(allMatrices,3)), 'FontAngle', 'italic');
    set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 11);

    colorBar = colorbar;
    colorBar.Label.String = 'Decoder accuracy';

end
