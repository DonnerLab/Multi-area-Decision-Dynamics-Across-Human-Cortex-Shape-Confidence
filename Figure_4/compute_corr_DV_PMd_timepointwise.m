function [All_corr_raw, All_corr_resid] = ...
    compute_corr_DV_PMd_timepointwise( ...
        subj_data, valid_areas, valid_ref, ...
        idx_t, choiceSign, sess)

%% ============================================================
% Correlation between PMd/M1 signal and area DVs
%
% OUTPUTS
%
% All_corr_raw
%   Correlation between PMd/M1 DV and raw area DV
%
% All_corr_resid
%   Correlation between PMd/M1 DV and PMd-projected area DV
%
% INPUTS
%
% subj_data
%   Structure containing decoder outputs
%
% valid_areas
%   Areas of interest (excluding PMd/M1 reference areas)
%
% valid_ref
%   Reference PMd/M1 areas used for projection:
%   typically {'6a','6d','4'}
%
% idx_t
%   Time window indices
%
% choiceSign
%   Trial-wise choice sign
%
% sess
%   Session labels
%
% NOTE
%
% DVs are first averaged across the selected time window.
% Then PMd/M1 is computed as the average across reference areas.
%
% We compute:
%
% 1. Raw correlation:
%    corr(PMd, DV)
%
% 2. Residual correlation:
%    corr(PMd, DV after projection-out of PMd)
%
%% ============================================================

nTrials = numel(choiceSign);
nAreas  = numel(valid_areas);

%% ============================================================
% COMPUTE WINDOW-AVERAGED DV FOR ALL AREAS
%% ============================================================

DV = nan(nTrials, nAreas);

for ai = 1:nAreas

    DV(:, ai) = mean( ...
        subj_data.(valid_areas{ai}).logodds(:, idx_t), ...
        2);

end

%% ============================================================
% COMPUTE PMd/M1 REFERENCE DV
%% ============================================================

ref = nan(nTrials, numel(valid_ref));

for ri = 1:numel(valid_ref)

    ref(:, ri) = mean( ...
        subj_data.(valid_ref{ri}).logodds(:, idx_t), ...
        2);

end

% PMd/M1 reference signal
DV_PMd = mean(ref, 2, 'omitnan');

%% ============================================================
% RAW CORRELATION
% corr(PMd, raw DV)
%% ============================================================

All_corr_raw = nan(nAreas, 1);

for ai = 1:nAreas

    keep = ...
        ~isnan(DV_PMd) & ...
        ~isnan(DV(:, ai));

    if sum(keep) > 2

        All_corr_raw(ai) = corr( ...
            DV_PMd(keep), ...
            DV(keep, ai), ...
            'rows', 'complete');

    end
end

%% ============================================================
% PROJECT OUT PMd/M1 (SESSION-WISE)
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

    % normalized PMd vector
    x1n = x1 / norm(x1);

    for ai = 1:nAreas

        % remove projection onto PMd axis
        DV_resid(idx, ai) = ...
            DV(idx, ai) ...
            - (x1n' * DV(idx, ai)) * x1n;

    end
end

%% ============================================================
% RESIDUAL CORRELATION
% corr(PMd, PMd-projected DV)
%% ============================================================

All_corr_resid = nan(nAreas, 1);

for ai = 1:nAreas

    keep = ...
        ~isnan(DV_PMd) & ...
        ~isnan(DV_resid(:, ai));

    if sum(keep) > 2

        All_corr_resid(ai) = corr( ...
            DV_PMd(keep), ...
            DV_resid(keep, ai), ...
            'rows', 'complete');

    end
end

end