%% =======================================================================
% Figure 2B Decoder Weights on Selected ROIs (Left + Right Hemisphere)
%=================================================

clear all
close all

%% =======================================================================
% Paths
%% =======================================================================

% Folder where this script is located
scriptFolder = fileparts(mfilename('fullpath'));

% Go up two folders 
projectFolder = fileparts(fileparts(scriptFolder));
%%
% Main folders
decoderFolder   = fullfile(projectFolder, 'data', 'brain_map_decoder');
decoderVXFolder = fullfile(projectFolder, 'data', 'decoder_output', 'all_VX');
functionsFolder = fullfile(projectFolder , 'functions');
giftiFolder     = fullfile(functionsFolder, 'gifti-1.6');

% Add required paths
addpath(decoderFolder)
addpath(decoderVXFolder)
addpath(functionsFolder)
addpath(giftiFolder)

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

%% =======================================================================
% Subject IDs
%% =======================================================================

subjectIDs = { ...
    'S01','S04','S06','S07','S12','S15','S17','S18', ...
    'S11','S05','S19','S20','S22','S10','S02','S08', ...
    'S24','S21','S23','S25'};

nSubjects = 20;
nSessions = 6;

% Selected ROIs from original code
selectedROIs = [27, 44, 47];

%% =======================================================================
% Load atlas
%% =======================================================================

load(fullfile(functionsFolder, 'Glasser_surface.mat'))

areaTable = readtable(fullfile(functionsFolder, 'Glasser_labels.csv'));
atlasLabels = areaTable.Label;

%% =======================================================================
% Plot setup
%% =======================================================================

figure

%% =======================================================================
% Loop across hemispheres
%% =======================================================================

for hemisphereSide = 1:2
    
    % Store vertex values for all subjects
    vertexValuesAllSubjects = nan(64984, nSubjects);
    
    %% ===================================================================
    % Loop across subjects
    %% ===================================================================
    
    for subject = 1:nSubjects
        
        subjectBetas = [];
        nVerticesPerSession = nan(1, nSessions);
        
        %% ===============================================================
        % First pass: determine minimum number of vertices across sessions
        %% ===============================================================
        
        for session = 1:nSessions
            
            fileName = [ ...
                char(subjectIDs(subject)), ...
                '_train_CHOICE_decoder_session', ...
                num2str(session), ...
                '_PCA_PMd_M1_SINGLETRIAL.mat'];
            
            load(fileName)
            
            nVerticesPerSession(session) = size(Prob_Performance, 3);
        end
        
        minVertices = min(nVerticesPerSession);
        
        %% ===============================================================
        % Second pass: load betas using same number of vertices
        %% ===============================================================
        
        for session = 1:nSessions
            
            fileName = [ ...
                char(subjectIDs(subject)), ...
                '_train_CHOICE_decoder_session', ...
                num2str(session), ...
                '_PCA_PMd_M1_SINGLETRIAL.mat'];
            
            load(fileName)
            
            subjectBetas(:,:,:,session) = ...
                Prob_Performance(:,:,1:minVertices);
        end
        
        meanBetas = nanmean(subjectBetas, 4);
        
        %% ===============================================================
        % Separate hemisphere
        %% ===============================================================
        
        if hemisphereSide == 1
            % Left hemisphere
            hemisphereBetas = mean( ...
                mean(meanBetas(1,561:601,1:n_vx_side(1)), 2), 1);
        else
            % Right hemisphere
            hemisphereBetas = mean( ...
                mean(meanBetas(1,561:601,...
                (n_vx_side(1)+1):n_vx_side(2)), 2), 1);
        end
        
        %% ===============================================================
        % Assign selected ROIs to cortical vertices
        %% ===============================================================
        
        subplot(1,2,hemisphereSide)
        
        for roi = selectedROIs
            
            if hemisphereSide == 1
                roiName = ['L_', char(allROIs(roi))];
            else
                roiName = ['R_', char(allROIs(roi))];
            end
            
            labelIndex = find(strcmp(atlasLabels, roiName));
            
            if ~isempty(labelIndex)
                
                colorMap = [brewermap(180, '*RdBu'); 1 1 1];
                
                atlasIndex = atlas.indexmax;
                vertexIndex = find(atlasIndex == labelIndex);
                
                % Keep same number of vertices across subjects
                sampledBetas = hemisphereBetas(:,:, ...
                    randi(size(hemisphereBetas,3), ...
                    [size(vertexIndex,1), 1]));
                
                for v = 1:size(vertexIndex,1)
                    vertexValuesAllSubjects(vertexIndex(v), subject) = ...
                        sampledBetas(:,:,v);
                end
            end
        end
    end
    
    %% ===================================================================
    % Average across subjects
    %% ===================================================================
    
    averageVertexValues = nanmean(vertexValuesAllSubjects, 2);
    
    %% ===================================================================
    % Hemisphere-specific plotting
    %% ===================================================================
    
    if hemisphereSide == 1
        hemi = 'L';
        surfaceType = 'inflated';
        plotValues = averageVertexValues(32493:end);
    else
        hemi = 'R';
        surfaceType = 'inflated';
        plotValues = averageVertexValues(1:32492);
    end
    
    %% ===================================================================
    % Load surface
    %% ===================================================================
    
    surfaceFile = fullfile( ...
        functionsFolder, ...
        ['S1200.' hemi '.' surfaceType '_MSMAll.32k_fs_LR.surf.gii']);
    
    surfaceMesh = gifti(surfaceFile);
    
    %% ===================================================================
    % Plot
    %% ===================================================================
    
    colormap(colorMap)
    
    % Make non-selected areas white/gray
    nanIndex = isnan(plotValues);
    plotValues(nanIndex) = 100;
    
    clim = [-0.0025 0.0025];
    
    trisurf( ...
        surfaceMesh.faces, ...
        surfaceMesh.vertices(:,1), ...
        surfaceMesh.vertices(:,2), ...
        surfaceMesh.vertices(:,3), ...
        plotValues, ...
        'edgealpha', 0)
    
    set(gca, 'clim', clim)
    set(gca, 'fontsize', 20)
    
    axis vis3d
    axis square
    axis off
    
    material dull
    lighting flat
    set(gcf, 'color', 'w')
    
    if hemisphereSide == 1
        view([270 0])
    else
        view([90 0])
    end
    
    xlim([-100 50])
    ylim([-105 70])
    zlim([-70 80])
    
    axis off
    camlight headlight
    
    hold on
end