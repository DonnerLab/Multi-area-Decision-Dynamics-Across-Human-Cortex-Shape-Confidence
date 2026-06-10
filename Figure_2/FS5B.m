%% =======================================================================
% V1 individual sample contrast decoder (Figure S5B)
%% =======================================================================

clear; close all;

% Paths
scriptFolder  = fileparts(mfilename('fullpath'));
projectFolder = fileparts(fileparts(scriptFolder));

dataFolder = fullfile(projectFolder, 'data', 'decoder_output', 'v1');

%% =======================================================================
% Subjects
%% =======================================================================

SubjectID = {'S01','S02','S04','S05','S07','S08','S11','S12', ...
             'S15','S17','S18','S10','S06','S19','S20', ...
             'S21','S22','S23','S24','S25'};

nSubjects = numel(SubjectID);
nSessions = 6;
nHemis    = 2;
nSamples  = 10;

%% =======================================================================
% Loop across hemispheres
%% =======================================================================

for hemi = 1:nHemis

    Decoder_kernel_all = [];

    for s = 1:nSubjects

        Decoder_ses = [];

        %% ----------------------------------------------------------------
        % Loop across sessions
        %% ----------------------------------------------------------------

        for ses = 1:nSessions

            if hemi == 1
                fileName = sprintf( ...
                    '%s_right_decoder_session%d_stimulus.mat', ...
                    SubjectID{s}, ses);
            else
                fileName = sprintf( ...
                    '%s_left_decoder_session%d_stimulus.mat', ...
                    SubjectID{s}, ses);
            end

            filePath = fullfile(dataFolder, fileName);

            if ~isfile(filePath)
                continue
            end

            load(filePath); % contains Decoder_Performance

            Decoder_ses = cat(3, Decoder_ses, Decoder_Performance);

        end

        if isempty(Decoder_ses)
            continue
        end

        %% ----------------------------------------------------------------
        % Average across sessions
        %% ----------------------------------------------------------------

        Decoder_mean = mean(Decoder_ses, 3, 'omitnan');
        Decoder_mean = Decoder_mean';

        Decoder_kernel_all(s,:,:) = Decoder_mean;

    end

    %% ===================================================================
    % Plot
    %% ===================================================================

    figure(hemi);
    hold on;

    colors = gray(11);

    for sam = 1:nSamples

        %% ---------------------------------------------------------------
        % Sample-specific plotting window
        %% ---------------------------------------------------------------

        sampleOnset = 201 + 40 * (sam - 1);

        if sam == 10
            sampleOffset = 301 + 40 * (sam - 1);
        else
            sampleOffset = 241 + 40 * (sam - 1);
        end

        plotWindow = (sampleOnset - 40):(sampleOffset + 40);

        %% ---------------------------------------------------------------
        % Plot average decoder precision
        %% ---------------------------------------------------------------

        y = mean( ...
            Decoder_kernel_all(:, plotWindow, sam), ...
            1, ...
            'omitnan');

        plot( ...
            plotWindow, ...
            y, ...
            'color', colors(sam,:), ...
            'linewidth', 3);

    end

    %% -------------------------------------------------------------------
    % Figure formatting
    %% -------------------------------------------------------------------

    ylabel('Decoder precision');
    xlabel('Time from stimulus onset');
    xticks([201 601]);
    xticklabels([0 1]);
    set(gca, 'fontsize', 10);

end