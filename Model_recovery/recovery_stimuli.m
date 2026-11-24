function [left, right, side, contrast] = recovery_stimuli(pool, cfg)
%RECOVERY_STIMULI  Synthetic stimuli generated like in the experiment (experiment.m + sample_contrast.m).
%   pool : vector of QUEST contrast differences to draw from (cfg.stim.contrast_mode = 'empirical'; ignored if 'fixed')
%   left, right : 10 x N contrast samples (times cfg.stim.scale);  side : +1 right stronger (correct), -1 left;
%   contrast : generative contrast difference used on each trial.
%
%   Per run (120 trials): side of the stronger grating = Markov chain with P(repeat) 0.5 (odd runs) or 0.8 (even runs);
%   10% of the trials are 'easy' (contrast difference 0.5).
%   Per trial: strong ~ N(baseline + c/2, sigma), weak ~ N(baseline - c/2, sigma), 10 samples each, values <0 -> 0.01,
%   values >1 -> 1, everything resampled until mean(strong) >= mean(weak).
sc  = cfg.stim;
ntr = sc.trials_per_run;
N   = sc.n_runs * ntr;
left = zeros(10, N);  right = zeros(10, N);  side = zeros(1, N);  contrast = zeros(1, N);
col = 0;
for r = 1:sc.n_runs
    if mod(r,2) == 1, Tm = [0.5 0.5; 0.5 0.5]; else, Tm = [0.8 0.2; 0.2 0.8]; end
    sd = zeros(1, ntr);
    sd(1) = randi(2);
    for t = 2:ntr
        sd(t) = find(cumsum(Tm(sd(t-1),:)) > rand(), 1);
    end
    sd = 2*sd - 3;                                   % 1 -> -1 (left stronger), 2 -> +1 (right stronger)
    easy = false(1, ntr);
    easy(randperm(ntr, round(ntr*sc.easy_fraction))) = true;
    for t = 1:ntr
        col = col + 1;
        if easy(t)
            c = sc.easy_contrast;
        elseif strcmp(sc.contrast_mode, 'empirical')
            c = pool(randi(numel(pool)));
        else
            c = sc.contrast_fixed;
        end
        c = min(1, max(0, c));
        [small, large] = local_sample_contrast(c, sc.sigma, sc.baseline);
        if sd(t) == -1
            left(:,col) = large';  right(:,col) = small';
        else
            left(:,col) = small';  right(:,col) = large';
        end
        side(col) = sd(t);  contrast(col) = c;
    end
end
left  = left  * sc.scale;
right = right * sc.scale;
end

function [small, large] = local_sample_contrast(contrast, sigma, baseline)
% identical to sample_contrast.m
large = [0, 0];  small = [0, 1];
while mean(large) < mean(small)
    large = randn(1,10)*sigma + contrast/2 + baseline;
    large(large < 0) = 0.01;  large(large > 1) = 1;
    small = randn(1,10)*sigma - contrast/2 + baseline;
    small(small < 0) = 0.01;  small(small > 1) = 1;
end
end
