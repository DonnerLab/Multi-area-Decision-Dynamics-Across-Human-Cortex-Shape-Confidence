function [xfit, info] = recovery_fit(f, dat, cfg)
%RECOVERY_FIT  Fit model f to dat.response/dat.conf with cost_function_model (same strategy as model_fit.m, cheaper):
%   best of cfg.n_starts random starts -> BADS with noisy-target handling.
%   Only the parameters used by the model are free; the others stay at cfg.nominal.
used = cfg.models(f).used;
nf   = numel(used);
lb = cfg.lb(used);  ub = cfg.ub(used);
lo = cfg.start_lo(used);  hi = cfg.start_hi(used);

expand = @(x) local_expand(x, used, cfg.nominal);
fun    = @(x) cost_function_model(f, expand(x), dat);

t0 = tic;

% 1) random starts
X0   = lo + (hi - lo).*rand(cfg.n_starts, nf);
nll0 = nan(cfg.n_starts, 1);
for i = 1:cfg.n_starts
    nll0(i) = fun(X0(i,:));
end
[~, ib] = min(nll0);
x0 = X0(ib,:);

% 2) noise level of the target at the starting point (for BADS)
rep = zeros(1, cfg.n_noise_evals);
for r = 1:cfg.n_noise_evals, rep(r) = fun(x0); end
noise = max(std(rep), 0.1);

% 3) plausible bounds around the start (as in model_fit.m, with a minimum width)
rg  = ub - lb;
w   = max(abs(0.6*x0), 0.05*rg);
plb = max(x0 - w, lb + 0.01*rg);
pub = min(x0 + w, ub - 0.01*rg);
x0  = min(max(x0, plb), pub);

options.Display             = 'off';
options.UncertaintyHandling = true;
options.NoiseSize           = noise;
if ~isempty(cfg.max_fun_evals)
    options.MaxFunEvals     = cfg.max_fun_evals;
end

try
    [xb, fval, exitflag] = bads(fun, x0, lb, ub, plb, pub, [], options);
catch ME
    warning('BADS failed (%s); returning best random start.', ME.message);
    xb = x0;  fval = NaN;  exitflag = -99;
end

xfit = expand(xb);
info = struct('x0',expand(x0), 'nll_start',min(nll0), 'noise',noise, ...
              'bads_fval',fval, 'exitflag',exitflag, 'fit_time',toc(t0));
end

function th = local_expand(x, used, nominal)
th = nominal;
th(used) = x;
end
