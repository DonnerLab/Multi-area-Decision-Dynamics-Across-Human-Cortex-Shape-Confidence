function [dat, truth] = recovery_generate(g, d, cfg)
%RECOVERY_GENERATE  Synthetic dataset d from generating model g (fully deterministic given g, d, cfg).
%   1. stimuli generated like in the experiment (recovery_stimuli)
%   2. true parameters: uniform in the start ranges + plausibility filter
%   3. exact probability of each (choice, confidence) category per trial from model_g_simulations
%   4. sample one outcome per trial -> training responses; sample a second, independent set on the same
%      stimuli -> held-out responses

nS = numel(cfg.subject_rows);
s  = cfg.subject_rows(mod(d-1, nS) + 1);

pool = [];
if strcmp(cfg.stim.contrast_mode, 'empirical')
    subj = recovery_subject(cfg, s);
    pool = abs(subj(:,23));  pool = pool(~isnan(pool));
    assert(~isempty(pool), 'Column 23 of All_Data is empty for subject %d: use cfg.stim.contrast_mode=''fixed''', s);
end

rng(10000*g + d, 'twister');
[left, right, side, contrast] = recovery_stimuli(pool, cfg);
dat.left = left;  dat.right = right;  dat.t = size(left,1);  dat.side = side;

% true parameters
ok = false;
for attempt = 1:cfg.max_resample
    theta = local_sample_theta(g, d, s, nS, cfg);
    Ps    = recovery_probs(g, theta, dat, 1);
    recovery_check_probs(Ps, g);
    ok = local_plausible(Ps, dat, cfg);
    if ok, break; end
end

% exact generative probabilities, then sample responses
P_true = recovery_probs(g, theta, dat, cfg.eval_reps);
recovery_check_probs(P_true, g);
[dat.response,      dat.conf]      = local_sample(P_true);
[dat.response_test, dat.conf_test] = local_sample(P_true);

truth = struct('theta',theta, 'P_true',P_true, 'subject',s, 'plausible',ok, 'attempts',attempt, 'contrast',contrast);
end

% ------------------------------------------------------------------------
function theta = local_sample_theta(g, d, s, nS, cfg)
% parameters used by model g drawn uniformly in the start ranges of model_fit.m
used  = cfg.models(g).used;
theta = cfg.nominal;
lo = cfg.start_lo(used);  hi = cfg.start_hi(used);
theta(used) = lo + (hi - lo).*rand(1, numel(used));
end

function ok = local_plausible(P, dat, cfg)
pr  = P(:,1) + P(:,2);                       % P(choose right)
ph  = P(:,1) + P(:,3);                       % P(high confidence)
cr  = double(dat.side(:) == 1);              % 1 if right is the correct side
acc = mean(pr.*cr + (1-pr).*(1-cr));
inr = @(v,r) (v >= r(1)) && (v <= r(2));
ok  = inr(acc, cfg.plaus.acc) && inr(mean(pr), cfg.plaus.pright) && inr(mean(ph), cfg.plaus.phigh);
end

function [resp, conf] = local_sample(P)
c    = cumsum(P, 2);
u    = rand(size(P,1),1) .* c(:,end);
k    = 1 + sum(u > c(:,1:3), 2);             % category 1..4
resp = double(k <= 2);                       % 1 = right
conf = 1 + mod(k,2);                         % categories 1,3 = high (2); 2,4 = low (1)
end
