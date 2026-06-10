function [Corr_matrix Corr_matrix_res] = compute_timepointwise_corr_matrix( ...
    subj_data, valid_areas, valid_ref, idx_t, choiceSign, sess)

%% ============================================================
% Compute correlation matrix across areas + PMd
%% ============================================================

nTrials = numel(choiceSign);
nAreas  = numel(valid_areas);
nT      = numel(idx_t);

Corr_all = nan(nAreas + 1, nAreas + 1, nT);

for tt = 1:nT

    t = idx_t(tt);

    %% --------------------------------------------------------
    % Area decision variables at single timepoint
    %% --------------------------------------------------------

    DV_mat = nan(nTrials, nAreas);

    for ai = 1:nAreas
        DV_mat(:, ai) = ...
            subj_data.(valid_areas{ai}).logodds(:, t);
    end

    %% --------------------------------------------------------
    % PMd reference (average across reference areas)
    %% --------------------------------------------------------

    ref_mat = nan(nTrials, numel(valid_ref));

    for ri = 1:numel(valid_ref)
        ref_mat(:, ri) = ...
            subj_data.(valid_ref{ri}).logodds(:, t);
    end

    DV_PMd = mean(ref_mat, 2, 'omitnan');

    %% --------------------------------------------------------
    % OPTIONAL: PMd projection removal
    %
    % Currently OFF to match your original code
    %% --------------------------------------------------------

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

        x1n = x1 / norm(x1);

        for ai = 1:nAreas
            DV_resid(idx, ai) = ...
                DV_mat(idx, ai) ...
                - (x1n' * DV_mat(idx, ai)) * x1n;
        end
    end

    %% --------------------------------------------------------
    % Correlation matrix for this timepoint
    %% --------------------------------------------------------

    Corr_all(:, :, tt) = ...
        corrcoef([DV_mat DV_PMd]);

    Corr_all_res(:, :, tt) = ...
        corrcoef([DV_resid DV_PMd]);
 

end

%% ------------------------------------------------------------
% Average across timepoints
%% ------------------------------------------------------------

Corr_matrix = mean(Corr_all, 3, 'omitnan');
Corr_matrix_res = mean(Corr_all_res, 3, 'omitnan');

end