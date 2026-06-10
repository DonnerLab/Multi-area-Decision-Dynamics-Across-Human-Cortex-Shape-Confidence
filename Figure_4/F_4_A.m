%% ============================================================
% Brain surface map of  decoder performance
%
%Figure 4A
%
%% ============================================================

clear;
close all;
clc;


%% ============================================================
% PATHS
%% ============================================================


scriptFolder  = fileparts(mfilename('fullpath'));
projectFolder = fileparts(fileparts(scriptFolder));

%% ---------------- DATA ----------------

data_path = fullfile( ...
    projectFolder, ...
    'data', ...
    'brain_map_decoder');

%% ---------------- FUNCTIONS ----------------

functions_path = fullfile( ...
    projectFolder, ...
     ...
    'functions');

gifti_path = fullfile( ...
    functions_path, ...
    'gifti-1.6');

addpath(data_path);
addpath(functions_path);
addpath(gifti_path);

%% ---------------- FILES ----------------

decoder_file_raw = fullfile( ...
    data_path, ...
    'all_areas_choice.mat');

decoder_file_proj = fullfile( ...
    data_path, ...
    'all_areas_choice_prj.mat');

surface_file = fullfile( ...
    functions_path, ...
    'S1200.L.inflated_MSMAll.32k_fs_LR.surf.gii');

atlas_file = fullfile( ...
    data_path, ...
    'Glasser_surface.mat');

labels_file = fullfile( ...
    data_path, ...
    'Glasser_labels.csv');

%% ============================================================
% ROI LIST (180 Glasser parcels)
%% ============================================================

All_ROI = {... 
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


nROI = numel(All_ROI);

%% ============================================================
% LOAD DECODER PERFORMANCE
%% ============================================================

Decoder_performance = importdata( ...
    fullfile(data_path, 'all_areas_choice.mat'));

% Average across sessions
Decoder_performance = nanmean(Decoder_performance, 2);

% Subject x ROI
Decoder_performance = reshape(Decoder_performance, 20, nROI);

% Mean decoder performance per ROI
Performance = nanmean(Decoder_performance, 1);

Perf_no_prj = Performance;

% Shift for permutation testing
Decoder_performance = Decoder_performance - min(Performance);

%% ============================================================
% PERMUTATION TEST + FDR
%% ============================================================

pvals = nan(1, nROI);

for i = 1:nROI

    pvals(i) = permutationTest( ...
        Decoder_performance(:, i), ...
        zeros(20,1), ...
        10000);

end

qvals = fdr_bh(pvals, 0.01);

% significant ROIs after FDR
significant = qvals;

%% ============================================================
% EXCLUDED REGIONS + BEST 25 ROIs
%% ============================================================

exclude_list = { ...
    '4','6a','6d','1','2','3a','3b'};

keep_mask     = ~ismember(All_ROI, exclude_list);
excluded_mask =  ismember(All_ROI, exclude_list);

% remove excluded areas from ranking
perf_tmp = Performance;
perf_tmp(~keep_mask) = -inf;

[~, idx_sorted] = sort(perf_tmp, 'descend');

best25_idx = idx_sorted(1:25);

best25_mask = false(1, nROI);
best25_mask(best25_idx) = true;

%% ============================================================
% LOAD GLASSER ATLAS
%% ============================================================

load('Glasser_surface.mat');

Areas = readtable('Glasser_labels.csv');
Areas = Areas.Label;


%% ============================================================
% COLORMAP
%% ============================================================

cmap_red = brewermap(256, 'Reds');

cmap = [ ...
    0.55 0.55 0.55;   % non-significant
    0.80 0.80 0.80;   % excluded
    cmap_red];

IDX_NONSIG  = 1;
IDX_EXCLUDE = 2;
IDX_OFFSET  = 2;

%% ============================================================
% NORMALIZE PERFORMANCE FOR COLOR MAPPING
%% ============================================================

perf_norm = Performance;


perf_norm = ...
    (perf_norm - 0.5) ...
    ./ (0.53 - min(Perf_no_prj));

%% ============================================================
% VERTEX COLOR ASSIGNMENT
%% ============================================================

% full cortex
VX = ones(64984,1) * IDX_NONSIG;

for i = 1:nROI

    ROIname = ['L_', All_ROI{i}];

    id = find(strcmp(Areas, ROIname));

    if isempty(id)
        continue
    end

    vidx = find(atlas.indexmax == id);

    if excluded_mask(i)

        VX(vidx) = IDX_EXCLUDE;

    elseif significant(i)

        idx = round(perf_norm(i) * 255) + 1 + IDX_OFFSET;
        VX(vidx) = idx;

    else

        VX(vidx) = IDX_NONSIG;

    end
end

%% ============================================================
% LEFT HEMISPHERE ONLY
%% ============================================================

VX_L = VX(32493:end);

%% ============================================================
% LOAD SURFACE
%% ============================================================

g = gifti(fullfile( ...
    functions_path, ...
    'S1200.L.inflated_MSMAll.32k_fs_LR.surf.gii'));

faces = g.faces;
verts = g.vertices;

%% ============================================================
% PLOT LEFT HEMISPHERE
%% ============================================================

figure('Color', 'w');

for view_i = 1:2

    subplot(2,1,view_i);

    colormap(cmap);

    trisurf( ...
        faces, ...
        verts(:,1), ...
        verts(:,2), ...
        verts(:,3), ...
        VX_L, ...
        'edgealpha', 0);

    hold on;

    caxis([1 size(cmap,1)]);

    %% --------------------------------------------------------
    % PARCEL BOUNDARIES
    %% --------------------------------------------------------

    for i = 1:nROI

        ROIname = ['L_', All_ROI{i}];

        id = find(strcmp(Areas, ROIname));

        if isempty(id)
            continue
        end

        vidx = find(atlas.indexmax == id);

        % left hemisphere indexing
        vidx = vidx(vidx >= 32493) - 32492;

        if isempty(vidx)
            continue
        end

        roi_mask = false(size(VX_L));
        roi_mask(vidx) = true;

        % faces fully inside ROI
        roi_faces = all(roi_mask(faces), 2);

        if ~any(roi_faces)
            continue
        end

        F = faces(roi_faces,:);

        % boundary edges
        E = [ ...
            F(:,[1 2]); ...
            F(:,[2 3]); ...
            F(:,[3 1])];

        E = sort(E, 2);

        [uE, ~, ic] = unique(E, 'rows');
        counts = accumarray(ic, 1);

        boundaryE = uE(counts == 1, :);

        %% --------------------------------------------
        % boundary color
        %% --------------------------------------------

        if best25_mask(i)

            col = [0 0 0];      % black
            lw  = 1.2;

        elseif excluded_mask(i)

            col = [0.3 0.3 0.3];
            lw  = 1.0;

        else

            continue

        end

        for e = 1:size(boundaryE,1)

            v = boundaryE(e,:);

            plot3( ...
                verts(v,1), ...
                verts(v,2), ...
                verts(v,3), ...
                'Color', col, ...
                'LineWidth', lw);

        end
    end

    %% --------------------------------------------------------
    % FIGURE STYLE
    %% --------------------------------------------------------

    axis off;
    shading interp;
    lighting gouraud;

    material([0.9 0.25 0 1 1]);

    camlight('headlight', 'infinite');
    camlight(90,0, 'infinite');
    camlight(-90,0, 'infinite');
    camlight(0,90, 'infinite');
    camlight(180,90, 'infinite');

    if view_i == 1
        view([90 0]);   % lateral
    else
        view([270 0]);  % medial
    end

    xlim([-100 50]);
    ylim([-105 70]);
    zlim([-70 80]);

end

