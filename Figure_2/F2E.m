%% =======================================================================
% Figure 2E
%
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

decoderFolder   = fullfile(projectFolder, 'data', 'decoder_output', 'all_VX');
functionsFolder = fullfile(projectFolder, 'functions');

addpath(functionsFolder)
addpath(fullfile(functionsFolder, 'gifti-1.6'))


%% =======================================================================
% SETTINGS
%% =======================================================================

subjects = { ...
    'S01','S04','S06','S07','S12','S15','S17','S18','S11','S05', ...
    'S19','S20','S22','S10','S02','S08','S24','S21','S23','S25'};

nSubjects = numel(subjects);
nSessions = 6;

% Atlas areas of PMd/M1
roiNames = {'4', '6d', '6a'};

% Last 100 ms of stimulus viewing (samples of the decoder output)
offsetWindow = 561:601;



% Glasser atlas (vertex labels)
load(fullfile(functionsFolder, 'Glasser_surface.mat'))
atlasTable  = readtable(fullfile(functionsFolder, 'Glasser_labels.csv'));
atlasLabels = atlasTable.Label;


%% =======================================================================
% LOOP OVER HEMISPHERES (1 = left, 2 = right)
%% =======================================================================

hemisphereNames = {'L', 'R'};

figure
colormap([brewermap(180, '*RdBu'); 1 1 1])

for hemisphere = 1:2

    vertexValues = nan(64984,nSubjects);

    for subject = 1:nSubjects

        %% ---------------------------------------------------------------
        % DECODER WEIGHTS OF THIS SUBJECT (average across sessions)
        %
        % Prob_Performance: 1 x time x vertices (decoder weights)
        % n_vx_side       : number of vertices of the left hemisphere and
        %                   of both hemispheres
        % Sessions are cut to the smallest number of vertices.
        %% ---------------------------------------------------------------

        sessionWeights = cell(1,nSessions);

        for session = 1:nSessions
            load(fullfile(decoderFolder, [subjects{subject} '_train_CHOICE_decoder_session' ...
                num2str(session) '_PCA_PMd_M1_SINGLETRIAL.mat']));
            sessionWeights{session} = Prob_Performance;
        end

        minVertices = min(cellfun(@(w) size(w,3), sessionWeights));

        subjectWeights = [];
        for session = 1:nSessions
            subjectWeights(:,:,:,session) = sessionWeights{session}(:,:,1:minVertices);
        end

        meanWeights = nanmean(subjectWeights,4);

        % Weights of this hemisphere, averaged over the last 100 ms
        if hemisphere == 1
            hemisphereVertices = 1:n_vx_side(1);
        else
            hemisphereVertices = (n_vx_side(1)+1):n_vx_side(2);
        end

        hemisphereWeights = squeeze(mean(meanWeights(1,offsetWindow,hemisphereVertices),2));


        %% ---------------------------------------------------------------
        % FILL THE ATLAS VERTICES OF THE 3 AREAS
        %% ---------------------------------------------------------------

        for roi = 1:numel(roiNames)

            labelIndex   = find(strcmp(atlasLabels, [hemisphereNames{hemisphere} '_' roiNames{roi}]));
            atlasVertices = find(atlas.indexmax == labelIndex);

            sampledVertices = randi(numel(hemisphereWeights), [numel(atlasVertices) 1]);

            vertexValues(atlasVertices,subject) = hemisphereWeights(sampledVertices);

        end

    end


    %% ===================================================================
    % FIGURE (Figure 2E) - AVERAGE ACROSS SUBJECTS
    %% ===================================================================

    averageValues = nanmean(vertexValues,2);

    if hemisphere == 1
        plotValues = averageValues(32493:end);
    else
        plotValues = averageValues(1:32492);
    end

    % Areas outside PMd/M1 in white
    plotValues(isnan(plotValues)) = 100;

    surfaceMesh = gifti(fullfile(functionsFolder, ...
        ['S1200.' hemisphereNames{hemisphere} '.inflated_MSMAll.32k_fs_LR.surf.gii']));

    subplot(1,2,hemisphere)

    trisurf( ...
        surfaceMesh.faces, ...
        surfaceMesh.vertices(:,1), ...
        surfaceMesh.vertices(:,2), ...
        surfaceMesh.vertices(:,3), ...
        plotValues, ...
        'edgealpha', 0)

    set(gca, 'clim', [-0.0025 0.0025])

    axis vis3d
    axis off
    if hemisphere == 1
        view([270 0])
    else
        view([90 0])
    end
    xlim([-100 50])
    ylim([-105 70])
    zlim([-70 80])

    set(gcf, 'color', 'w')
    material dull
    lighting flat
    camlight headlight

end
