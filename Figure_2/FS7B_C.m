%% =======================================================================
% Figure S7B-C
%
% Decoding of the individual contrast samples from V1
%
% S7B: left-hemifield contrast samples decoded from right V1
% S7C: right-hemifield contrast samples decoded from left V1
%
% Decoding precision = cross-validated correlation between decoded and
% presented contrast
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

% V1 decoder output
decoderFolder = fullfile(projectFolder, 'data', 'decoder_output', 'v1');


%% =======================================================================
% SETTINGS
%% =======================================================================

subjects = { ...
    'S01','S02','S04','S05','S07','S08','S11','S12','S15','S17', ...
    'S18','S10','S06','S19','S20','S21','S22','S23','S24','S25'};

nSubjects = numel(subjects);
nSessions = 6;
nSamples  = 10;

% Decoder of each panel: right V1 (left contrasts), left V1 (right contrasts)
hemisphereNames = {'right', 'left'};
panelTitles     = {'Decoding of left contrast', 'Decoding of right contrast'};

% One gray level per sample position (sample 1 = black)
sampleColors = gray(11);


%% =======================================================================
% LOOP OVER HEMISPHERES
%% =======================================================================

for hemisphere = 1:2

    % subjects x time x sample position
    decodingPrecision = [];

    for subject = 1:nSubjects

        % Decoder_Performance: sample position x time; average across sessions
        sessionPrecision = [];

        for session = 1:nSessions

            fileName = fullfile(decoderFolder, sprintf('%s_%s_decoder_session%d_stimulus.mat', ...
                subjects{subject}, hemisphereNames{hemisphere}, session));

            if isfile(fileName)
                load(fileName);
                sessionPrecision = cat(3, sessionPrecision, Decoder_Performance);
            end

        end

        decodingPrecision(subject,:,:) = nanmean(sessionPrecision,3)';

    end


    %% ===================================================================
    % FIGURE (Figure S7B / S7C)
    %
    % Time base of the decoder output: 400 samples per second, stimulus
    % onset at sample 201, one contrast sample every 40 samples (100 ms);
    % each sample position is shown from 100 ms before to 200 ms after
    % its onset (350 ms for the last sample)
    %% ===================================================================

    figure(hemisphere)
    hold on

    for sample = 1:nSamples

        sampleOnset = 201 + 40*(sample - 1);

        if sample == nSamples
            sampleOffset = sampleOnset + 100;
        else
            sampleOffset = sampleOnset + 40;
        end

        plotWindow = (sampleOnset - 40):(sampleOffset + 40);

        plot(plotWindow, nanmean(decodingPrecision(:,plotWindow,sample),1), ...
            'color', sampleColors(sample,:), 'linewidth', 3);

    end

    xticks([201 601])
    xticklabels([0 1])
    xlabel('Time from stimulus onset (s)')
    ylabel({'Decoding precision'; '(correlation coeff.)'})
    title(panelTitles{hemisphere}, 'FontWeight', 'normal')
    set(gca, 'fontsize', 15)

end
