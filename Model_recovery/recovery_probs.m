function P = recovery_probs(m, theta, dat, reps)
%RECOVERY_PROBS  Probability of the 4 (choice, confidence) categories per trial under model m.
%   Columns: [R/high, R/low, L/high, L/low] = the 'Confidence' output of model_m_simulations(pars, data).
%   reps > 1 averages independent calls (equivalent to 1000*reps simulations per trial; the models
%   normalise confidence per simulation across trials, so averaging calls is exact).
if nargin < 4, reps = 1; end
d.left  = dat.left;
d.right = dat.right;
d.t     = size(dat.left, 1);
fh = str2func(sprintf('model_%d_simulations', m));
P = 0;
for r = 1:reps
    [~, Pr, ~] = fh(theta, d);
    Pr = double(Pr);
    if size(Pr,2) ~= 4, Pr = Pr'; end
    P = P + Pr;
end
P = P / reps;
end
