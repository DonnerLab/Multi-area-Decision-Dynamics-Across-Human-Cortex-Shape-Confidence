%% =======================================================================
% Figure 2A Choice Decoder Performance Across Cortical Areas
%
%=========================================

clear all
close all

%% =======================================================================
% Paths
%% =======================================================================

% Folder where this script is located
scriptFolder = fileparts(mfilename('fullpath'));

% Go up two folders 
projectFolder = fileparts(fileparts(scriptFolder));

% Data folders
decoderFolder = fullfile(projectFolder, 'data', 'brain_map_decoder');
functionsFolder = fullfile(projectFolder, 'functions');
giftiFolder = fullfile(functionsFolder, 'gifti-1.6');

% Add required paths
addpath(decoderFolder)
addpath(functionsFolder)
addpath(giftiFolder)

%% =======================================================================
% Load decoder performance
%% =======================================================================

decoderData = importdata(fullfile(decoderFolder, 'all_areas_choice.mat'));

%% =======================================================================
% ROI labels
%% =======================================================================

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

nSubjects = 20;
nAreas = 180;

%% =======================================================================
% Average decoder performance
%% =======================================================================

% Original code:
% mean across repetitions, then subtract chance level (0.5)

decoderPerformance = nanmean(decoderData, 2) - 0.5;
decoderPerformance = reshape(decoderPerformance, nSubjects, nAreas);

% Mean performance across subjects for each area
meanPerformance = nan(1, nAreas);

for area = 1:nAreas
    meanPerformance(area) = nanmean(decoderPerformance(:, area));
end

%% =======================================================================
% Load Glasser atlas
%% =======================================================================

load(fullfile(functionsFolder, 'Glasser_surface.mat'))

areaTable = readtable(fullfile(functionsFolder, 'Glasser_labels.csv'));
atlasLabels = areaTable.Label;

%% =======================================================================
% Assign decoder values to cortical vertices
%% =======================================================================

vertexValues = nan(64984, 1);

for area = 1:nAreas
    
    roiName = ['L_', char(allROIs(area))];
    
    labelIndex = find(strcmp(atlasLabels, roiName));
    
    if ~isempty(labelIndex)
        
        atlasIndex = atlas.indexmax;
        vertexIndex = find(atlasIndex == labelIndex);
        
        vertexValues(vertexIndex) = meanPerformance(area);
        
    end
end

%% =======================================================================
% Colormap
%% =======================================================================

colorMap = [brewermap(180, '*RdBu'); 0.5 0.5 0.5];

%% =======================================================================
% Left hemisphere only
%% =======================================================================

% Original code uses only left hemisphere vertices
leftHemisphereValues = vertexValues(32493:end);

maxValue = max(abs(leftHemisphereValues));

% Replace NaNs with gray color index
nanIndex = isnan(leftHemisphereValues);
leftHemisphereValues(nanIndex) = maxValue + 0.5;

%% =======================================================================
% Load cortical surface
%% =======================================================================

hemisphere = 'L';
surfaceType = 'inflated';

surfaceFile = fullfile( ...
    functionsFolder, ...
    ['S1200.' hemisphere '.' surfaceType '_MSMAll.32k_fs_LR.surf.gii'] ...
);

surfaceMesh = gifti(surfaceFile);

%% =======================================================================
% Plot brain surface
%% =======================================================================

figure
colormap(colorMap)

clim = [-0.03 0.03];

trisurf( ...
    surfaceMesh.faces, ...
    surfaceMesh.vertices(:,1), ...
    surfaceMesh.vertices(:,2), ...
    surfaceMesh.vertices(:,3), ...
    leftHemisphereValues, ...
    'edgealpha', 0 ...
)

set(gca, 'clim', clim)
set(gca, 'fontsize', 20)

axis vis3d
axis square
axis off

material dull
lighting flat
set(gcf, 'color', 'w')

view([90 0])

xlim([-100 50])
ylim([-105 70])
zlim([-70 80])

%% =======================================================================
% Final rendering settings
%% =======================================================================

material([0.6 0.4 0 5 1])

shading interp
lighting gouraud

% Symmetric lighting to reduce shadows
camlight('headlight')
camlight('left')
camlight('right')