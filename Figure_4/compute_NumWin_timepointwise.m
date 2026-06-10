function NumWin = compute_NumWin_timepointwise( ...
    subj_data, valid_areas, valid_ref, idx_t, choiceSign, sess)

% ============================================================
% Compute number of winning areas after PMd projection removal
%
% Winning area:
% sign(area DV residual) == subject choice
%
% Steps:
% 1. Average DV across selected time window
% 2. Compute PMd reference DV
% 3. Regress PMd projection out 
% 4. Count winning areas
% ============================================================

nTrials = numel(choiceSign);
nAreas  = numel(valid_areas);

%% ============================================================
% AREA DVs (average across selected time window)
%% ============================================================

DV_mat = nan(nTrials, nAreas);

for ai = 1:nAreas
    DV_mat(:, ai) = mean( ...
        subj_data.(valid_areas{ai}).logodds(:, idx_t), ...
        2 ...
    );
end

%% ============================================================
% PMd REFERENCE DV
%% ============================================================

ref_mat = nan(nTrials, numel(valid_ref));

for ri = 1:numel(valid_ref)
    ref_mat(:, ri) = mean( ...
        subj_data.(valid_ref{ri}).logodds(:, idx_t), ...
        2 ...
    );
end

DV_PMd = mean(ref_mat, 2, 'omitnan');

%% ============================================================
% REMOVE PMd PROJECTION (session-wise)
%% ============================================================

DV_resid = DV_mat;

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
        DV_resid(idx, ai) = DV_mat(idx, ai) ...
            - (x1n' * DV_mat(idx, ai)) * x1n;

    end
end

%% ============================================================
% COUNT WINNING AREAS
%% ============================================================

win_mask = sign(DV_resid) == choiceSign;

% remove exact zeros
win_mask(DV_resid == 0) = false;

NumWin = sum(win_mask, 2);

end