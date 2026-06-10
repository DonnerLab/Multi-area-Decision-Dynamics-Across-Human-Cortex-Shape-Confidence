function [AvgWinDV, AvgLoseDV] = ...
    compute_AvgWin_AvgLose_timepointwise( ...
    subj_data, valid_areas, valid_ref, idx_t, choiceSign, sess)

% ============================================================
% Compute average winning DV and average losing DV
% after removing PMd projection
%
% Winning area:
%     sign(area DV residual) == choice
%
% Losing area:
%     sign(area DV residual) ~= choice
%
% Steps:
% 1. Average DV across selected time window
% 2. Compute PMd reference DV
% 3. Remove PMd projection 
% 4. Identify winning vs losing areas
% 5. Compute mean absolute residual DV for:
%       - winning areas
%       - losing areas
% ============================================================

nTrials = numel(choiceSign);
nAreas  = numel(valid_areas);

%% ============================================================
% AREA DVs (average across selected time window)
%% ============================================================

DV = nan(nTrials, nAreas);

for ai = 1:nAreas

    DV(:, ai) = mean( ...
        subj_data.(valid_areas{ai}).logodds(:, idx_t), ...
        2 ...
    );

end

%% ============================================================
% PMd REFERENCE DV
%% ============================================================

ref = nan(nTrials, numel(valid_ref));

for ri = 1:numel(valid_ref)

    ref(:, ri) = mean( ...
        subj_data.(valid_ref{ri}).logodds(:, idx_t), ...
        2 ...
    );

end

DV_PMd = mean(ref, 2, 'omitnan');

%% ============================================================
% REMOVE PMd PROJECTION (session-wise)
%% ============================================================

DV_resid = DV;

for ss = unique(sess)'

    idx = sess == ss;

    if sum(idx) < 3
        continue
    end

    x1 = DV_PMd(idx);

    if norm(x1) == 0
        continue
    end

    % normalize PMd vector
    x1n = x1 / norm(x1);

    for ai = 1:nAreas

        % remove projection onto PMd axis
        DV_resid(idx, ai) = DV(idx, ai) ...
            - (x1n' * DV(idx, ai)) * x1n;

    end
end

%% ============================================================
% DEFINE WINNING AND LOSING AREAS
%% ============================================================

% winning = same sign as choice
win_mask = sign(DV_resid) == choiceSign;

% losing = opposite sign and not zero
lose_mask = sign(DV_resid) ~= choiceSign & DV_resid ~= 0;

%% ============================================================
% COMPUTE MEAN ABSOLUTE RESIDUAL DV
%% ============================================================

% winning DV magnitude
W = abs(DV_resid);
W(~win_mask) = NaN;

% losing DV magnitude
L = abs(DV_resid);
L(~lose_mask) = NaN;

AvgWinDV = mean(W, 2, 'omitnan');
AvgLoseDV = mean(L, 2, 'omitnan');

end