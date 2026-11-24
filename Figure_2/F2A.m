%% =======================================================================
% Figure 2A
%
% Whole-cortex maps of choice and confidence decoding accuracy at the end
% of stimulus viewing (left hemisphere, 180 areas of the Glasser atlas)
%
% Decoder accuracy of each area: data/brain_map_decoder/
%   all_areas_choice.mat, all_areas_confidence.mat
% (all_areas_confidence.mat: functions/decoder/decode_confidence_bilateral.m,
% then functions/decoder/make_brain_map_matrix.m)
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

% Data and functions (atlas, surfaces, gifti toolbox)
decoderFolder   = fullfile(projectFolder, 'data', 'brain_map_decoder');
functionsFolder = fullfile(projectFolder, 'functions');

addpath(functionsFolder)
addpath(fullfile(functionsFolder, 'gifti-1.6'))


%% =======================================================================
% SETTINGS
%% =======================================================================

nSubjects = 20;
nAreas    = 180;

% Area names, in the order of the decoder output
allROIs = { ...
    'V1','V2','V3','V4','V6','V3A','V7','IPS1','V3B','V6A','V8','FFC','PIT','VMV1','VMV3','VMV2','VVC',...
    'MST','LO1','LO2','MT','PH','V4t','FST','V3CD','LO3','4','3b','1','2','3a','5m','5mv',...
    '23c','5L','24dd','24dv','SCEF','6ma','6mp','FEF','PEF','55b','6d','6v','6r','6a',...
    '43','OP4','OP1','OP2-3','FOP1','A1','52','RI','PFcm','PBelt','MBelt','LBelt','TA2','STGa',...
    'A5','STSda','STSdp','STSvp','A4','STSva','PoI2','FOP4','MI','Pir','AVI','AAIC','FOP3','FOP2','PoI1',...
    'Ig','FOP5','PI','EC','PreS','H','PeEc','PHA1','PHA3','TF','PHA2','TGd','TE1a','TE1p','TE2a','TE2p',...
    'PHT','TGv','TE1m','PSL','STV','TPOJ1','TPOJ2','TPOJ3','7Pm','7AL','7Am','7PL','7PC','LIPv','VIP',...
    'MIP','LIPd','AIP','PFt','PGp','IP2','IP1','IP0','PFop','PF','PFm','PGi','PGs','RSC','POS2','PCV','7m',...
    'POS1','23d','v23ab','d23ab','31pv','ProS','DVT','31pd','31a','p24pr','33pr','a24pr','p32pr','a24',...
    'd32','8BM','p32','10r','9m','10v','25','s32','pOFC','a32pr','p24','47m','10d','a10p','10pp','11l','13l',...
    'OFC','47s','p10p','44','45','47l','a47r','IFJa','IFJp','IFSp','IFSa','p47r','SFL','8Av','8Ad','8BL',...
    '9p','8C','p9_46v','46','a9_46v','9_46d','9a','i6_8','s6_8'};

% Glasser atlas (vertex labels) and inflated left-hemisphere surface
load(fullfile(functionsFolder, 'Glasser_surface.mat'))
atlasTable  = readtable(fullfile(functionsFolder, 'Glasser_labels.csv'));
atlasLabels = atlasTable.Label;

surfaceMesh = gifti(fullfile(functionsFolder, 'S1200.L.inflated_MSMAll.32k_fs_LR.surf.gii'));


%% =======================================================================
% ONE MAP FOR EACH DECODER (choice, confidence)
%% =======================================================================

decoderNames = {'choice', 'confidence'};

for decoder = 1:2

    %% ===================================================================
    % DECODER ACCURACY OF EACH AREA (group average, minus chance)
    %% ===================================================================

    decoderData = importdata(fullfile(decoderFolder, ['all_areas_' decoderNames{decoder} '.mat']));

    % Mean across columns (sessions) minus chance level, subjects x areas
    decoderAccuracy = reshape(nanmean(decoderData,2) - 0.5, nSubjects, nAreas);

    meanAccuracy = nanmean(decoderAccuracy,1);


    %% ===================================================================
    % ASSIGN THE ACCURACY OF EACH AREA TO ITS VERTICES
    %% ===================================================================

    vertexValues = nan(64984,1);

    for area = 1:nAreas

        labelIndex = find(strcmp(atlasLabels, ['L_' allROIs{area}]));

        if ~isempty(labelIndex)
            vertexValues(atlas.indexmax == labelIndex) = meanAccuracy(area);
        end

    end

    % Left-hemisphere vertices; vertices outside the atlas in gray
    leftHemisphereValues = vertexValues(32493:end);
    leftHemisphereValues(isnan(leftHemisphereValues)) = max(abs(leftHemisphereValues)) + 0.5;


    %% ===================================================================
    % FIGURE (Figure 2A)
    %% ===================================================================

    figure(decoder)
    colormap([brewermap(180, '*RdBu'); 0.5 0.5 0.5])

    trisurf( ...
        surfaceMesh.faces, ...
        surfaceMesh.vertices(:,1), ...
        surfaceMesh.vertices(:,2), ...
        surfaceMesh.vertices(:,3), ...
        leftHemisphereValues, ...
        'edgealpha', 0)

    set(gca, 'clim', [-0.03 0.03])
    title(['Decoding of ' decoderNames{decoder}])

    axis vis3d
    axis off
    view([90 0])
    xlim([-100 50])
    ylim([-105 70])
    zlim([-70 80])

    set(gcf, 'color', 'w')
    material([0.6 0.4 0 5 1])
    shading interp
    lighting gouraud
    camlight('headlight')
    camlight('left')
    camlight('right')


end
