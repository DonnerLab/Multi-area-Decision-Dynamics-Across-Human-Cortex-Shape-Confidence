function recovery_job(task_id)
%% =======================================================================
% One cluster task = one fit: generating model g, dataset d, fitted
% model f (task list: recovery_tasklist.m)
%
% Saves data/model/recovery/rec_gGG_dDD_fFF.mat with the true and fitted
% parameters, in-sample log-likelihood, held-out log-likelihood
% (independent responses to the same stimuli), AIC and BIC
%% =======================================================================

cfg = recovery_setup();
addpath(cfg.code_dir);  addpath(cfg.bads_path);  addpath(cfg.models_dir);
if ~exist(cfg.out_dir,'dir'), mkdir(cfg.out_dir); end

T = recovery_tasklist(cfg);
if task_id < 1 || task_id > size(T,1)
    error('task_id %d out of range (1..%d)', task_id, size(T,1));
end
g = T(task_id,1);  d = T(task_id,2);  f = T(task_id,3);
for mm = unique([g f])
    assert(exist(sprintf('model_%d_simulations', mm), 'file') == 2, ...
           'model_%d_simulations.m not found: set cfg.models_dir in recovery_setup.m', mm);
end

outfile = fullfile(cfg.out_dir, sprintf('rec_g%02d_d%02d_f%02d.mat', g, d, f));
if exist(outfile,'file')
    return
end

t0 = tic;
[dat, truth] = recovery_generate(g, d, cfg);

% keep the synthetic dataset (same for all fitted models of this dataset)
if ~exist(cfg.data_dir,'dir'), mkdir(cfg.data_dir); end
dfile = fullfile(cfg.data_dir, sprintf('data_g%02d_d%03d.mat', g, d));
if ~exist(dfile,'file')
    dtmp = [tempname(cfg.data_dir) '.mat'];
    save(dtmp, 'dat', 'truth', '-v7');
    movefile(dtmp, dfile);
end

% fail fast if the fitted model cannot produce valid probabilities
dsmall.left = dat.left(:,1:20);  dsmall.right = dat.right(:,1:20);
recovery_check_probs(recovery_probs(f, cfg.nominal, dsmall, 1), f);

rng(100000*f + 1000*g + d, 'twister');
[xfit, fi] = recovery_fit(f, dat, cfg);

% final, low-noise evaluation of the fitted parameters
P_eval   = recovery_probs(f, xfit, dat, cfg.eval_reps);
LL_train = recovery_loglik(P_eval, dat.response,      dat.conf);
LL_test  = recovery_loglik(P_eval, dat.response_test, dat.conf_test);
k        = numel(cfg.models(f).used);
n        = numel(dat.response);

R = struct();
R.g = g;  R.d = d;  R.f = f;  R.k = k;  R.n_trials = n;
R.theta_true = truth.theta;  R.theta_fit = xfit;
R.used_true  = cfg.models(g).used;  R.used_fit = cfg.models(f).used;
R.LL_train = LL_train;  R.LL_test = LL_test;
R.AIC = -2*LL_train + 2*k;
R.BIC = -2*LL_train + k*log(n);
R.LL_oracle_train = recovery_loglik(truth.P_true, dat.response,      dat.conf);
R.LL_oracle_test  = recovery_loglik(truth.P_true, dat.response_test, dat.conf_test);
R.plausible = truth.plausible;  R.subject = truth.subject;
R.fitinfo = fi;  R.total_time = toc(t0);

tmp = fullfile(cfg.out_dir, sprintf('tmp_model_%d.mat', task_id));
save(tmp, '-struct', 'R', '-v7');
movefile(tmp, outfile);                   % no half-written result files
end
