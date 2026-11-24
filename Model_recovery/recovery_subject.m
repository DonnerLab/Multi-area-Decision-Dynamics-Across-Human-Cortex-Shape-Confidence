function subj = recovery_subject(cfg, s)
%RECOVERY_SUBJECT  All trials of subject s (all drug conditions), without trials lacking a confidence report.
%   Same aggregation as model_fit.m. Columns: 1 response, 2:11 left contrasts, 12:21 right contrasts,
%   22 confidence-1, 23 generative mean, 24 RT, 25:26 session id.
S = load(cfg.data_file, 'All_Data');
subj = [];
for dd = 1:3
    subj = [subj; cell2mat(S.All_Data(s,dd))]; %#ok<AGROW>
end
subj(isnan(subj(:,22)), :) = [];
end
