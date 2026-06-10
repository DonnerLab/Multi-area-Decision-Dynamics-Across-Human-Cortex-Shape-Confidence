function [AvgWinDV, AvgLoseDV, ...
          AvgWinDV_raw, AvgLoseDV_raw, ...
          AvgWinDV_nabs, AvgLoseDV_nabs, ...
          AvgWinDV_raw_nabs, AvgLoseDV_raw_nabs, ...
          N_win_areas] = ...
    compute_AvgWin_AvgLose_timepointwise_correlations( ...
        subj_data, valid_areas, valid_ref, ...
        idx_t, choiceSign, sess)

%% ============================================================
% Compute Avg Winning DV + Avg Losing DV
%
% OUTPUTS
%
% AvgWinDV            : projected winning DV (absolute value)
% AvgLoseDV           : projected losing DV (absolute value)
%
% AvgWinDV_raw        : raw winning DV (absolute value)
% AvgLoseDV_raw       : raw losing DV (absolute value)
%
% AvgWinDV_nabs       : projected winning DV (signed, no abs)
% AvgLoseDV_nabs      : projected losing DV (signed, no abs)
%
% AvgWinDV_raw_nabs   : raw winning DV (signed, no abs)
% AvgLoseDV_raw_nabs  : raw losing DV (signed, no abs)
%
% N_win_areas         : number of winning areas per trial
%% ============================================================

nTrials = numel(choiceSign);
nAreas  = numel(valid_areas);
nT      = numel(idx_t);

%% ============================================================
% STORAGE
%% ============================================================

Win_tp           = nan(nTrials, nT);
Lose_tp          = nan(nTrials, 1);

Win_tp_raw       = nan(nTrials, 1);
Lose_tp_raw      = nan(nTrials, 1);

Win_tp_nabs      = nan(nTrials, 1);
Lose_tp_nabs     = nan(nTrials, 1);

Win_tp_raw_nabs  = nan(nTrials, 1);
Lose_tp_raw_nabs = nan(nTrials, 1);

N_win_areas      = nan(nTrials, 1);

%% ============================================================
% MAIN LOOP
%% ============================================================

for tt = 1
    % Keeping tt = 1 exactly as in the original code

    %% --------------------------------------------------------
    % Area DVs (window-averaged)
    %% --------------------------------------------------------

    DV = nan(nTrials, nAreas);

    for ai = 1:nAreas
        DV(:, ai) = mean( ...
            subj_data.(valid_areas{ai}).logodds(:, idx_t), ...
            2);
    end

    %% --------------------------------------------------------
    % PMd reference (window-averaged)
    %% --------------------------------------------------------

    ref = nan(nTrials, numel(valid_ref));

    for ri = 1:numel(valid_ref)
        ref(:, ri) = mean( ...
            subj_data.(valid_ref{ri}).logodds(:, idx_t), ...
            2);
    end

    DV_PMd = mean(ref, 2, 'omitnan');

    %% --------------------------------------------------------
    % Remove PMd/M1 projection (session-wise)
    %% --------------------------------------------------------

    DV_resid = DV;
    DV_raw   = DV;

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
                DV(idx, ai) ...
                - (x1n' * DV(idx, ai)) * x1n;

        end
    end

    %% --------------------------------------------------------
    % PROJECTED WIN / LOSE (ABSOLUTE VALUE)
    %% --------------------------------------------------------

    win  = sign(DV_resid) == choiceSign;
    lose = sign(DV_resid) ~= choiceSign & DV_resid ~= 0;

    W = abs(DV_resid);
    L = abs(DV_resid);

    W(~win)  = NaN;
    L(~lose) = NaN;

    N_win_areas(:, tt) = sum(win, 2);

    %% --------------------------------------------------------
    % RAW WIN / LOSE (ABSOLUTE VALUE)
    %% --------------------------------------------------------

    win  = sign(DV_raw) == choiceSign;
    lose = sign(DV_raw) ~= choiceSign & DV_raw ~= 0;

    W_raw = abs(DV_raw);
    L_raw = abs(DV_raw);

    W_raw(~win)  = NaN;
    L_raw(~lose) = NaN;

    %% --------------------------------------------------------
    % PROJECTED WIN / LOSE (SIGNED, NO ABS)
    %% --------------------------------------------------------

    win  = sign(DV_resid) == choiceSign;
    lose = sign(DV_resid) ~= choiceSign & DV_resid ~= 0;

    W_nabs = DV_resid;
    L_nabs = DV_resid;

    W_nabs(~win)  = NaN;
    L_nabs(~lose) = NaN;

    %% --------------------------------------------------------
    % RAW WIN / LOSE (SIGNED, NO ABS)
    %% --------------------------------------------------------

    win  = sign(DV_raw) == choiceSign;
    lose = sign(DV_raw) ~= choiceSign & DV_raw ~= 0;

    W_raw_nabs = DV_raw;
    L_raw_nabs = DV_raw;

    W_raw_nabs(~win)  = NaN;
    L_raw_nabs(~lose) = NaN;

    %% --------------------------------------------------------
    % Trial-wise averages
    %% --------------------------------------------------------

    Win_tp(:, tt)  = mean(W, 2, 'omitnan');
    Lose_tp(:, tt) = mean(L, 2, 'omitnan');

    Win_tp_raw(:, tt)  = mean(W_raw, 2, 'omitnan');
    Lose_tp_raw(:, tt) = mean(L_raw, 2, 'omitnan');

    Win_tp_nabs(:, tt)  = mean(W_nabs, 2, 'omitnan');
    Lose_tp_nabs(:, tt) = mean(L_nabs, 2, 'omitnan');

    Win_tp_raw_nabs(:, tt)  = mean(W_raw_nabs, 2, 'omitnan');
    Lose_tp_raw_nabs(:, tt) = mean(L_raw_nabs, 2, 'omitnan');

end

%% ============================================================
% FINAL AVERAGES ACROSS TIME
%% ============================================================

AvgWinDV  = mean(Win_tp, 2, 'omitnan');
AvgLoseDV = mean(Lose_tp, 2, 'omitnan');

AvgWinDV_raw  = mean(Win_tp_raw, 2, 'omitnan');
AvgLoseDV_raw = mean(Lose_tp_raw, 2, 'omitnan');

AvgWinDV_nabs  = mean(Win_tp_nabs, 2, 'omitnan');
AvgLoseDV_nabs = mean(Lose_tp_nabs, 2, 'omitnan');

AvgWinDV_raw_nabs  = mean(Win_tp_raw_nabs, 2, 'omitnan');
AvgLoseDV_raw_nabs = mean(Lose_tp_raw_nabs, 2, 'omitnan');

end